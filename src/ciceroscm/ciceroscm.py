"""
CICEROSCM
"""

import logging
import os

import numpy as np
import pandas as pd

from ._utils import cut_and_check_pamset, write_results_tsv
from .component_factory_functions import create_thermal_model
from .concentrations_emissions_handler import ConcentrationsEmissionsHandler
from .input_handler import InputHandler
from .make_plots import plot_output1
from .pub_utils import get_first_key

LOGGER = logging.getLogger(__name__)


FORC_OUTPUT_LIST = [
    "Total_forcing",
    "Solar_forcing",
    "Volcanic_forcing_NH",
    "Volcanic_forcing_SH",
]


class CICEROSCM:
    """
    Main ciceroscm class

    Takes care of and routes calls to diffusion
    model and concentrations emissions handler
    Takes results, and should be the module
    interacted with from other programs

    Attributes
    ----------
    cfg : dict
            Normalized configuration dictionary with the
            simulation year range, optional forcing-run
            settings, and the selected thermal and carbon
            cycle model names.
        thermal_model_class : type
            Thermal model class selected from the current
            configuration and instantiated for each run.
        rf : np.ndarray or pd.DataFrame
            External forcing input used for forcing-only runs.
            Present only when ``cfg["rf_run"]`` is true.
    ce_handler : ciceroscm.ConcentrationsEmissionsHandler
                   Concentrations and emissions handler for
                   concentration- or emissions-driven runs.
                   Present only when ``cfg["rf_run"]`` is false.
        feedback_list : list[str] or None
            Feedback variables retrieved from the concentrations
            and emissions handler. Present only when
            ``cfg["rf_run"]`` is false.
    results : dict
                Results dictionary for thermal-model output.
                When ``results_as_dict`` is requested, handler
                output such as emissions, concentrations, and
                component-wise forcings is merged into it as well.
        _sun_arr, _volc_n_arr, _volc_s_arr, _volc_n_mean,
        _volc_s_mean, _rf_luc_arr : np.ndarray
            Preloaded annual forcing inputs used during the model
            year loop. These arrays back solar, volcanic, and land-use
            forcing contributions.
    """

    # pylint: disable=too-many-instance-attributes
    # Consider whether this can be cut back later

    def __init__(self, cfg):
        """
        Intialise CICEROSCM

        Starting by picking out the part of cfg that
        are needed as the class cfgs.
        Then we either read in forc_file for forcing
        run, or check file availability and set start,
        end and emissions start for concentrations or
        emissions run.
        Then make empty dictionary for results, read in
        solar and volcanic forcing and initialise other
        output arrays

        Parameters
        ----------
        cfg : dict
              Configurations containing inputs about class
              such as a forcing file for forcing run,
              locations of files to use for concentration
              or emission runs, and start and end of run etc.

        """
        self.cfg = cut_and_check_pamset(
            {"nystart": 1750, "nyend": 2100, "emstart": 1850, "idtm": 24},
            cfg,
            {"carbon_cycle_model": "default", "thermal_model": "default"},
        )
        cfg.update(self.cfg)

        input_handler = InputHandler(cfg)
        self.cfg["rf_run"] = input_handler.optional_pam("forc")
        self.cfg["thermal_model"] = input_handler.thermal_model(self.cfg)
        self.cfg["carbon_cycle_model"] = input_handler.carbon_model(self.cfg)

        self.thermal_model_class = create_thermal_model(self.cfg["thermal_model"])

        #        print("Thermal Model=" + self.cfg["thermal_model"])

        if self.cfg["rf_run"]:
            self.rf = input_handler.get_data("forc")
        else:
            # cfg = check_inputfiles(cfg)
            pamset_emiconc = {}

            pamset_emiconc["carbon_cycle_model"] = self.cfg["carbon_cycle_model"]
            pamset_emiconc["emstart"] = self.cfg["emstart"]
            pamset_emiconc["nystart"] = self.cfg["nystart"]
            pamset_emiconc["nyend"] = self.cfg["nyend"]
            pamset_emiconc["idtm"] = self.cfg["idtm"]

            for key, value in cfg.items():
                if key in ["rs_function", "rb_function"]:
                    pamset_emiconc[key] = value
            self.ce_handler = ConcentrationsEmissionsHandler(
                input_handler, pamset_emiconc
            )
            self.feedback_list = self.ce_handler.get_feedback_list()
        self.results = {}
        # Solar / volcanic / LUC forcing as numpy arrays for the year loop;
        # these inputs are static for the run and the InputHandler returns
        # them as numpy directly, indexed positionally by (yr - nystart).
        # Sun and LUC are 1D; volcanic data is 2D (rows x hemispheric
        # subgrid) so we keep the per-year row for the thermal model and
        # pre-average for the annual output.
        self._sun_arr = input_handler.get_data("rf_sun")
        self._volc_n_arr = input_handler.get_data("rf_volc_n")
        self._volc_s_arr = input_handler.get_data("rf_volc_s")
        self._volc_n_mean = self._volc_n_arr.mean(axis=1)
        self._volc_s_mean = self._volc_s_arr.mean(axis=1)
        self._rf_luc_arr = input_handler.get_data("rf_luc")
        self.initialise_output_arrays()

    def initialise_output_arrays(self):
        """
        Initialise dict with arrays to hold data for run

        Dictionary for all results from upwelling diffusion
        model outputs is initialised with empty arrays
        """
        output_variables = (
            list(self.thermal_model_class.get_output_dict_thermal().keys())
            + FORC_OUTPUT_LIST
        )
        for output in output_variables:
            self.results[output] = np.zeros(self.cfg["nyend"] - self.cfg["nystart"] + 1)

    def forc_set(self, yr):
        """
        Read the forcing for this year

        Getting a year, getting the forcung for this year and
        adding solar forcing.

        Parameters
        ----------
        yr : int
          Year for which to read out data

        Returns
        -------
        tuple
            ``(fn, fs, forc, w_aero)``. The first three entries are the
            hemispheric and total forcing for the year (with solar
            added). ``w_aero`` is the magnitude-weighted aerosol
            forcing fraction used by the pattern-effect machinery; it
            comes from the optional ``w_aero`` column on the forcing
            DataFrame and defaults to ``0.0`` when the column is
            absent (graceful fallback for forcing files that do not
            carry per-agent decomposition).
        """
        row_index = yr - self.cfg["nystart"]
        # Add support for other forcing formats
        if isinstance(self.rf, np.ndarray):
            # Add luc albedo later
            forc = self.rf[row_index]  # + self.rf_luc.iloc[row_index][0]
            fn = forc
            fs = forc
            w_aero = 0.0
        else:
            forc = self.rf["total"][yr]
            fn = self.rf["FORC_NH"][yr]
            fs = self.rf["FORC_SH"][yr]
            w_aero = (
                float(self.rf["w_aero"][yr]) if "w_aero" in self.rf.columns else 0.0
            )
        forc = forc + self._sun_arr[row_index]
        fn = fn + self._sun_arr[row_index]
        fs = fs + self._sun_arr[row_index]
        return fn, fs, forc, w_aero

    def add_year_data_to_output(self, values, forc, index):
        """
        Add single year output to output arrays

        Add all the outputs from a single year run
        of upwelling diffusion model to output arrays

        Parameters
        ----------
        values : dict
              Output from upwelling diffusion model
        forc : float
            Total forcing for this year
        index : int
             Index equalling year number in the possible years
        """
        for output, name in self.thermal_model_class.get_output_dict_thermal().items():
            self.results[output][index] = values[name]
        self.results["Total_forcing"][index] = forc
        self.results["Solar_forcing"][index] = self._sun_arr[index]
        self.results["Volcanic_forcing_NH"][index] = self._volc_n_mean[index]
        self.results["Volcanic_forcing_SH"][index] = self._volc_s_mean[index]

    def _run(
        self,
        cfg,
        pamset_udm=None,
        pamset_emiconc=None,
        pamset_carbon=None,
        make_plot=False,
    ):  # pylint: disable=too-many-arguments, too-many-positional-arguments, too-many-locals, too-many-branches
        """
        Run CICEROSCM

        Setting off a full model run. Starting by
        intialising output arrays, and udm_model and
        resetting ConcentrationEmissionsHandler for a new run
        Then looping over year by year converting emissions
        and concnetrations to forcings if applicable
        and then running the upwelling diffusion model
        Finally writing results to file

        Parameters
        ----------
        cfg : dict
           Dictionary with run specific configurations
        pamset_udm : dict
                  Parameter set for udm model
        pamset_emiconc : dict
                      Parameter set for concentrations
                      emissions handler
        pamset_carbon : dict
            Parameter set for carbon cycle module
        make_plot : bool
            Whether to output plots automatically
        """
        self.initialise_output_arrays()
        # Setting up thermal model with parameters
        # udm = UpwellingDiffusionModel(pamset_udm)

        udm = self.thermal_model_class(pamset_udm)

        # Pattern-mediated feedback: each thermal model owns the
        # lambda_eff = lambda_0 + delta_lambda_aero * w_aero formula + delta_lambda_pdo
        # internally; the driver only forwards w_aero each year

        values = None
        if not self.cfg["rf_run"]:
            self.ce_handler.reset_with_new_pams(pamset_emiconc, pamset_carbon)
        for yr in range(self.cfg["nystart"], self.cfg["nyend"] + 1):
            if not self.cfg["rf_run"]:
                if values is not None and self.feedback_list is not None:
                    self.ce_handler.emi2conc(
                        yr, {key: values.get(key, 0) for key in self.feedback_list}
                    )
                else:
                    self.ce_handler.emi2conc(yr)
                forc, fn, fs, w_aero = self.ce_handler.conc2forc(
                    yr,
                    self._rf_luc_arr[yr - self.cfg["nystart"]],
                    self._sun_arr[yr - self.cfg["nystart"]],
                )

            else:
                fn, fs, forc, w_aero = self.forc_set(yr)

            values = udm.energy_budget(
                fn,
                fs,
                self._volc_n_arr[yr - self.cfg["nystart"]],
                self._volc_s_arr[yr - self.cfg["nystart"]],
                w_aero,
                yr - self.cfg["nystart"],
            )
            self.add_year_data_to_output(values, forc, yr - self.cfg["nystart"])

        if make_plot:
            plot_output1(cfg, self.results, self.cfg["nystart"], self.cfg["nyend"])
        if ("results_as_dict" in cfg) and cfg["results_as_dict"]:
            if not self.cfg["rf_run"]:
                self.results.update(
                    self.ce_handler.add_results_to_dict(
                        cfg,
                        feedback_dict_series={
                            key: self.results[
                                get_first_key(udm.get_output_dict_thermal(), key)
                            ]
                            for key in self.feedback_list
                        },
                    )
                )
        else:
            if not self.cfg["rf_run"]:
                self.ce_handler.write_output_to_files(
                    cfg,
                    feedback_dict_series={
                        key: self.results[
                            get_first_key(udm.get_output_dict_thermal(), key)
                        ]
                        for key in self.feedback_list
                    },
                    make_plot=make_plot,
                )

            self.write_data_to_file(cfg)

    def write_data_to_file(self, pamset):
        """
        Write results to files after run

        Writing results from upwelling diffusion model to file
        Formatting and organising in ocean heat content (ohc),
        radiative imbalance (rib), temperature related (temp),
        and forcing (forc) files are as in original fortran
        implementation. Forcing file is only outputted here
        if the run is a forcing run. Otherwise the forcing results
        writing is handled by the ConcentrationsEmissionsHandler

        Parameters
        ----------
        pamset : dict
              parameterset with details on where to write results

        """
        if "output_folder" in pamset:
            # Make os independent?
            outdir = os.path.join(os.getcwd(), pamset["output_folder"])
        else:
            outdir = os.path.join(os.getcwd(), "output")

        indices = np.arange(self.cfg["nystart"], self.cfg["nyend"] + 1)
        df_ohc = pd.DataFrame(
            data={
                "Year": indices,
                "OHC700": self.results["OHC700"],
                "OHCTOT": self.results["OHCTOT"],
            }
        )
        list_rib = ["RIB_glob", "RIB_N", "RIB_S"]
        df_rib = pd.DataFrame(data={"Year": indices})
        for vari in list_rib:
            df_rib[vari] = self.results[vari]
        list_temp = [
            "dT_glob",
            "dT_NH",
            "dT_SH",
            "dT_glob_air",
            "dT_NH_air",
            "dT_SH_air",
            "dT_glob_sea",
            "dT_NH_sea",
            "dT_SHsea",
        ]
        df_temp = pd.DataFrame(data={"Year": indices})
        if "output_prefix" in pamset:
            filename_start = pamset["output_prefix"]
        else:
            filename_start = "output"
        for vari in list_temp:
            df_temp[vari] = self.results[vari]

        write_results_tsv(df_ohc, os.path.join(outdir, f"{filename_start}_ohc.txt"))
        write_results_tsv(df_rib, os.path.join(outdir, f"{filename_start}_rib.txt"))
        write_results_tsv(df_temp, os.path.join(outdir, f"{filename_start}_temp.txt"))
        df_sunvolc = pd.DataFrame(
            data={
                "Year": indices,
                "Solar_forcing": self.results["Solar_forcing"],
                "Volcanic_forcing_NH": self.results["Volcanic_forcing_NH"],
                "Volcanic_forcing_SH": self.results["Volcanic_forcing_SH"],
            }
        )
        write_results_tsv(
            df_sunvolc, os.path.join(outdir, f"{filename_start}_sunvolc.txt")
        )
        if self.cfg["rf_run"]:
            df_forc = pd.DataFrame(
                data={"Year": indices, "Total_forcing": self.results["Total_forcing"]}
            )
            write_results_tsv(
                df_forc, os.path.join(outdir, f"{filename_start}_forc.txt")
            )
