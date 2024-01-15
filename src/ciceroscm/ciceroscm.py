"""
CICEROSCM
"""

import logging
import os

import numpy as np
import pandas as pd

from ._utils import cut_and_check_pamset
from .concentrations_emissions_handler import ConcentrationsEmissionsHandler
from .input_handler import InputHandler
from .make_plots import plot_output1
from .upwelling_diffusion_model import UpwellingDiffusionModel

LOGGER = logging.getLogger(__name__)


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
          Configurations dictionary with startyear,
          endyear, emissions start plus optional
          configs, like a forcing file path for pure
          forcing runs, else input files needed for
          concentrations and emission handler,
          a parameter to include volcanic and solar
          forcing, one to indicate a pure concentrations
          run, perturbation files etc.
    ce_handler : ciceroscm.ConcentrationsEmissionsHandler
                     Concentrations emissions handler for the
                     model. It can be reset for multiple runs
                     but needs to be used with the same year range
                     and with either Concentrations or emissions
                     run consistently
    results : dict
                  Results dictionary which stores results from
                  upwelling diffusion model. If a results_as_dict
                  is sent as True to the run call, the results
                  from the Concentrations emissions handler
                  (emissions, componentwise forcings and concentrations)
                  will also be added to this dictionary and can be
                  accessed from there
    rf_volc_sun : dict
                  With solar and volcanic forcing read from standard files
                  or set to zero according to configurations
    rf_luc : pd.Dataframe
             Dataframe of land use albedo forcing with years as index
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
            {"nystart": 1750, "nyend": 2100, "emstart": 1850, "idtm": 24}, cfg
        )
        cfg.update(self.cfg)
        input_handler = InputHandler(cfg)
        self.cfg["rf_run"] = input_handler.optional_pam("forc")
        if self.cfg["rf_run"]:
            self.rf = input_handler.get_data("forc")
        else:
            # cfg = check_inputfiles(cfg)
            pamset_emiconc = {}
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
        self.results = {}
        # Reading in solar and volcanic forcing
        self.rf_volc_sun = {
            "volc_n": input_handler.get_data("rf_volc_n"),
            "volc_s": input_handler.get_data("rf_volc_s"),
            "sun": input_handler.get_data("rf_sun"),
        }

        # Add support for sending filename in cfg
        self.rf_luc = input_handler.get_data("rf_luc")
        self.initialise_output_arrays()

    def initialise_output_arrays(self):
        """
        Initialise dict with arrays to hold data for run

        Dictionary for all results from upwelling diffusion
        model outputs is initialised with empty arrays
        """
        output_variables = [
            "OHC700",
            "OHCTOT",
            "RIB_glob",
            "RIB_N",
            "RIB_S",
            "dT_glob",
            "dT_NH",
            "dT_SH",
            "dT_glob_air",
            "dT_NH_air",
            "dT_SH_air",
            "dT_glob_sea",
            "dT_NH_sea",
            "dT_SHsea",
            "Total_forcing",
            "Solar_forcing",
            "Volcanic_forcing_NH",
            "Volcanic_forcing_SH",
        ]
        for output in output_variables:
            self.results[output] = np.zeros(self.cfg["nyend"] - self.cfg["nystart"] + 1)

    def forc_set(self, yr, rf_sun):
        """
        Read the forcing for this year

        Getting a year, getting the forcung for this year and
        adding solar forcing.

        Parameters
        ----------
        yr : int
          Year for which to read out data
        rf_sun : pandas.Dataframe
              Dataframe with solar forcing to  add to the other
              forcings.

        Returns
        -------
        float
             The total forcing for the year, including solar
             forcing is added.
        """
        row_index = yr - self.cfg["nystart"]
        # Add support for other forcing formats
        if isinstance(self.rf, np.ndarray):
            # Add luc albedo later
            forc = self.rf[row_index]  # + self.rf_luc.iloc[row_index][0]
            fn = forc
            fs = forc
        else:
            forc = self.rf["total"][yr]
            fn = self.rf["FORC_NH"][yr]
            fs = self.rf["FORC_SH"][yr]
        forc = forc + rf_sun.iloc[row_index, 0]
        fn = fn + rf_sun.iloc[row_index, 0]
        fs = fs + rf_sun.iloc[row_index, 0]
        return fn, fs, forc

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
        simple_outputs = ["OHC700", "OHCTOT"]
        for output in simple_outputs:
            self.results[output][index] = values[output]
        outputs_dict = {
            "RIB_glob": "RIB",
            "RIB_N": "RIBN",
            "RIB_S": "RIBS",
            "dT_glob": "dtemp",
            "dT_NH": "dtempnh",
            "dT_SH": "dtempsh",
            "dT_glob_air": "dtemp_air",
            "dT_NH_air": "dtempnh_air",
            "dT_SH_air": "dtempsh_air",
            "dT_glob_sea": "dtemp_sea",
            "dT_NH_sea": "dtempnh_sea",
            "dT_SHsea": "dtempsh_sea",
        }
        for output, name in outputs_dict.items():
            self.results[output][index] = values[name]
        self.results["Total_forcing"][index] = forc
        self.results["Solar_forcing"][index] = self.rf_volc_sun["sun"].iloc[index, 0]
        self.results["Volcanic_forcing_NH"][index] = np.mean(
            np.array(self.rf_volc_sun["volc_n"].iloc[index, :])
        )
        self.results["Volcanic_forcing_SH"][index] = np.mean(
            np.array(self.rf_volc_sun["volc_s"].iloc[index, :])
        )

    def _run(
        self, cfg, pamset_udm={}, pamset_emiconc={}, make_plot=False
    ):  # pylint: disable=dangerous-default-value
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
        """
        self.initialise_output_arrays()
        # Setting up UDM
        udm = UpwellingDiffusionModel(pamset_udm)
        if not self.cfg["rf_run"]:
            pamset_emiconc["emstart"] = self.cfg["emstart"]
            pamset_emiconc["nystart"] = self.cfg["nystart"]
            pamset_emiconc["nyend"] = self.cfg["nyend"]
            self.ce_handler.reset_with_new_pams(pamset_emiconc)
        for yr in range(self.cfg["nystart"], self.cfg["nyend"] + 1):
            if not self.cfg["rf_run"]:
                self.ce_handler.emi2conc(yr)
                forc, fn, fs = self.ce_handler.conc2forc(
                    yr,
                    self.rf_luc.iloc[yr - self.cfg["nystart"], 0],
                    self.rf_volc_sun["sun"].iloc[yr - self.cfg["nystart"], 0],
                )

            else:
                fn, fs, forc = self.forc_set(yr, self.rf_volc_sun["sun"])
            values = udm.energy_budget(
                fn,
                fs,
                np.array(self.rf_volc_sun["volc_n"].iloc[yr - self.cfg["nystart"], :]),
                np.array(self.rf_volc_sun["volc_s"].iloc[yr - self.cfg["nystart"], :]),
            )
            self.add_year_data_to_output(values, forc, yr - self.cfg["nystart"])

        if make_plot:
            plot_output1(cfg, self.results, self.cfg["nystart"], self.cfg["nyend"])
        if ("results_as_dict" in cfg) and cfg["results_as_dict"]:
            if not self.cfg["rf_run"]:
                self.results.update(self.ce_handler.add_results_to_dict())
        else:
            if not self.cfg["rf_run"]:
                self.ce_handler.write_output_to_files(cfg, make_plot)

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

        df_ohc.to_csv(
            os.path.join(outdir, f"{filename_start}_ohc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_rib.to_csv(
            os.path.join(outdir, f"{filename_start}_rib.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_temp.to_csv(
            os.path.join(outdir, f"{filename_start}_temp.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_sunvolc = pd.DataFrame(
            data={
                "Year": indices,
                "Solar_forcing": self.results["Solar_forcing"],
                "Volcanic_forcing_NH": self.results["Volcanic_forcing_NH"],
                "Volcanic_forcing_SH": self.results["Volcanic_forcing_SH"],
            }
        )
        df_sunvolc.to_csv(
            os.path.join(outdir, f"{filename_start}_sunvolc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        if self.cfg["rf_run"]:
            df_forc = pd.DataFrame(
                data={"Year": indices, "Total_forcing": self.results["Total_forcing"]}
            )
            df_forc.to_csv(
                os.path.join(outdir, f"{filename_start}_forc.txt"),
                sep="\t",
                index=False,
                float_format="%.5e",
            )
