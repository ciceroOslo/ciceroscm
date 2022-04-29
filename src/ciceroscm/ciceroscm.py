"""
CICEROSCM
"""
import logging
import os

import numpy as np
import pandas as pd

from ._utils import cut_and_check_pamset
from .concentrations_emissions_handler import ConcentrationsEmissionsHandler
from .upwelling_diffusion_model import UpwellingDiffusionModel

LOGGER = logging.getLogger(__name__)

default_data_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "default_data"
)


def check_inputfiles(cfg):
    """
    Check whether input files are present or not

    Checking configuration dictionary to see whether
    it necessary files for concentrations or
    emissions run are present.
    If natural emissions files are not found in cfg,
    standard file location is used.

    Parameters
    ----------
    cfg : dict
       Configurations dictionary which should contain
       locations of necessary files.

    Returns
    -------
    dict
        cfg possible augmented with standard locations
        for natural emissions files

    Raises
    ------
    FileNotFoundError
         If files are not found
    """
    if not os.path.exists(cfg["gaspamfile"]):
        raise FileNotFoundError(
            f"Concentration input file {cfg['concentrations_file']} not found"
        )
    if not os.path.exists(cfg["concentrations_file"]):
        raise FileNotFoundError(
            f"Concentration input file {cfg['concentrations_file']} not found"
        )
    if not os.path.exists(cfg["emissions_file"]):
        raise FileNotFoundError(
            f"Emission input file {cfg['emissions_file']} not found"
        )
    if "nat_ch4_file" not in cfg:
        LOGGER.warning(
            "Did not find prescribed nat_ch4_file or none was given. Looking in standard path",
        )
        cfg["nat_ch4_file"] = os.path.join(
            os.getcwd(), "input_OTHER", "NATEMIS", "natemis_ch4.txt"
        )
    if not os.path.exists(cfg["nat_ch4_file"]):
        raise FileNotFoundError(
            f"Natural emission input file {cfg['nat_ch4_file']} not found"
        )
    if "nat_n2o_file" not in cfg:
        LOGGER.warning(
            "Did not find prescribed nat_n2o_file or none was given. Looking in standard path",
        )
        cfg["nat_n2o_file"] = os.path.join(
            os.getcwd(), "input_OTHER", "NATEMIS", "natemis_n2o.txt"
        )
    if not os.path.exists(cfg["nat_n2o_file"]):
        raise FileNotFoundError(
            f"Natural emission input file {cfg['nat_n2o_file']} not found"
        )
    return cfg


def read_forc(forc_file):
    """
    Read in forcing from forc_file

    Read in forcing file to dataframe, couple of options
    depending on file formatting

    Parameters
    ----------
    forc_file : str
             Full path of forcing file to be read

    Returns
    -------
    ndarray
           Forcing data, or possibly a  pandas.Dataframe if
           data is organized in several components

    """
    components = False
    with open(forc_file, "r", encoding="utf8") as fread:
        first_line = fread.readline()
        if first_line[:4].lower() == "year":
            components = True
    if not components:
        df_forc = np.loadtxt(forc_file)
    else:
        # Decide on formatting for this
        df_forc = pd.read_csv(forc_file, index_col=0)
    return df_forc


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

        Raises
        ------
        FileNotFoundError
            If forcing file is not found when forcing run is chosen
        """
        self.cfg = cut_and_check_pamset(
            {"nystart": 1750, "nyend": 2100, "emstart": 1850}, cfg
        )
        rf_run = False
        if "forc_file" in cfg:
            rf_run = True
            if not os.path.exists(cfg["forc_file"]):
                raise FileNotFoundError(
                    f"Forcing input file {cfg['forc_file']} not found"
                )
            self.rf = read_forc(cfg["forc_file"])
        else:
            cfg = check_inputfiles(cfg)
            pamset_emiconc = {}
            pamset_emiconc["emstart"] = self.cfg["emstart"]
            pamset_emiconc["nystart"] = self.cfg["nystart"]
            pamset_emiconc["nyend"] = self.cfg["nyend"]
            self.ce_handler = ConcentrationsEmissionsHandler(cfg, pamset_emiconc)

        self.cfg["rf_run"] = rf_run
        self.results = {}
        # Reading in solar and volcanic forcing
        self.rf_volc_sun = self.read_in_volc_and_sun(cfg)

        # Add support for sending filename in cfg
        self.rf_luc = self.read_data_on_year_row(
            os.path.join(default_data_dir, "IPCC_LUCalbedo.txt")
        )
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
            "dSL(m)",
            "dSL_thermal(m)",
            "dSL_ice(m)",
            "Total_forcing",
        ]
        for output in output_variables:
            self.results[output] = np.zeros(self.cfg["nyend"] - self.cfg["nystart"] + 1)

    def read_data_on_year_row(self, volc_datafile):
        """
        Read in data from file with no headers


        Read in data from file with no headers where
        each year is a row. Typically this is the format for
        volcano and solar data. The years are taken to be
        the years from the defined startyear and endyear

        Parameters
        ----------
        volc_datafile : str
                     Path of file to be read

        Returns
        -------
        pandas.Dataframe
                        Dataframe containing the data with the years as
                        indices
        """
        indices = np.arange(self.cfg["nystart"], self.cfg["nyend"] + 1)
        nrows = len(indices)
        if self.cfg["nystart"] > 1750:
            skiprows = self.cfg["nystart"] - 1750
            df_data = pd.read_csv(
                volc_datafile,
                header=None,
                skiprows=skiprows,
                nrows=nrows,
                delim_whitespace=True,
            )
        else:
            df_data = pd.read_csv(
                volc_datafile, header=None, nrows=nrows, delim_whitespace=True
            )

        df_data.set_axis(labels=indices, inplace=True)
        return df_data

    def read_in_volc_and_sun(self, cfg):
        """
        Read in solar and volcanic forcing and return them

        Read in solar or volcanic forcing if this is chosen
        otherwise produce empty dataframes that can be used
        instead. If solar and volcanic forcing is added, a
        hemispherically dependent addition is added to the
        volcanic part, to adjust for lack of spin up.

        Parameters
        ----------
        cfg : dict
           Dictionary containing configurations on whether to use
           solar and volcanic forcing or not.

        Returns
        -------
        dict
            Containing the dataframes for hemispheric volcanic
            forcings and solar forcing
        """
        if "sunvolc" in cfg and cfg["sunvolc"] == 1:
            # Possibly change to allow for other files
            # And for SH to differ from NH
            rf_volc_n = self.read_data_on_year_row(
                os.path.join(default_data_dir, "meanVOLCmnd_ipcc_NH.txt")
            )
            rf_volc_s = rf_volc_n
            # Test, adjust for not including spin up. See Gregory et al.
            # Se regneark.
            rf_volc_n = rf_volc_n + 0.371457071
            rf_volc_s = rf_volc_s + 0.353195076
            rf_sun = self.read_data_on_year_row(
                os.path.join(default_data_dir, "solar_IPCC.txt")
            )
        # Add support for sending filename in cfg
        else:
            indices = np.arange(self.cfg["nystart"], self.cfg["nyend"] + 1)
            rf_volc_n = pd.DataFrame(
                data=np.zeros((self.cfg["nyend"] - self.cfg["nystart"] + 1, 12)),
                index=indices,
                columns=range(12),
            )
            rf_volc_s = rf_volc_n
            rf_sun = pd.DataFrame(
                data={0: np.zeros(self.cfg["nyend"] - self.cfg["nystart"] + 1)}
            )
        return {"volc_n": rf_volc_n, "volc_s": rf_volc_s, "sun": rf_sun}

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
        else:
            forc = self.rf["total"][row_index]
        forc = forc + rf_sun.iloc[row_index, 0]
        return forc

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
        self.results["dSL(m)"][index] = values["deltsl"][0] + values["deltsl"][1]
        self.results["dSL_ice(m)"][index] = values["deltsl"][1]
        self.results["dSL_thermal(m)"][index] = values["deltsl"][0]
        self.results["Total_forcing"][index] = forc

    def _run(
        self, cfg, pamset_udm={}, pamset_emiconc={}
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
                    self.rf_luc.loc[yr, 0],
                    self.rf_volc_sun["sun"].loc[yr - self.cfg["nystart"], 0],
                )

            else:
                forc = self.forc_set(yr, self.rf_volc_sun["sun"])
                fs = forc
                fn = forc
            values = udm.energy_budget(
                fn,
                fs,
                self.rf_volc_sun["volc_n"].iloc[yr - self.cfg["nystart"], :],
                self.rf_volc_sun["volc_s"].iloc[yr - self.cfg["nystart"], :],
            )
            self.add_year_data_to_output(values, forc, yr - self.cfg["nystart"])

        if "results_as_dict" in cfg and cfg["results_as_dict"]:
            self.results.update(self.ce_handler.add_results_to_dict())
        else:
            if not self.cfg["rf_run"]:
                self.ce_handler.write_output_to_files(cfg)

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
            "dSL(m)",
            "dSL_thermal(m)",
            "dSL_ice(m)",
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
