"""
Module to take care of reading in inputs
wither from files or as direct objects
"""

import logging
import os

import numpy as np
import pandas as pd

# from ._utils import check_numeric_pamset
from ._utils import cut_and_check_pamset

LOGGER = logging.getLogger(__name__)

default_data_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "default_data"
)


def read_components(filename):
    """
    Read in components to be considered

    Read the gas_pam files and rename a few headers to
    make it easier to use

    Parameters
    ----------
    filename : int
            path to gaspamfile

    Returns
    -------
    pandas.Dataframe
        Dataframe with gases, and various info on them
    """
    df_gas = pd.read_csv(filename, delim_whitespace=True, index_col=0)
    df_gas.rename(
        columns={"TAU1(YEARS)": "TAU1", "NATURAL_EMISSIONS": "NAT_EM"}, inplace=True
    )
    return df_gas


def read_natural_emissions(filename, component, startyear=1750, endyear=2500):
    """
    Read in single column natural emissions file

    A natural emissions file with data on a single column is read in
    to a pandas Dataframe. A an index of corresponding years is
    generated and added. Data in input is assumed to be yearly.

    Parameters
    ----------
    filename : str
               path to file with emissions
    component : str
                Name of component the emissions are for
    startyear: int
               Startyear for emissions
    endyear : int
              Endyear for emissions

    Returns
    -------
    pd.Dataframe
                Dataframe of natural emissions for component with years as index
    """
    df_natemis = pd.read_csv(filename, header=None, names=[component], index_col=False)
    df_natemis["year"] = np.arange(startyear, endyear + 1)
    df_natemis = df_natemis.set_index("year")
    return df_natemis


def read_inputfile(input_file, cut_years=False, year_start=1750, year_end=2100):
    """
    Read input from emissions or concentrations file

    Parameters
    ----------
    input_file : str
              Path to file to be read in
    cut_years : bool
             If unused years are to be cut, this option should
             be set to True, default is False
    year_start : int
             Start year for relevant input, default is 1750
    year_end : int
             End year for relevant input, default is 2100

    Returns
    -------
    pandas.Dataframe
        Dataframe with the intput from the file, possibly cut
        to relevant years
    """
    df_input = pd.read_csv(
        input_file, delim_whitespace=True, index_col=0, skiprows=[1, 2, 3]
    )
    if cut_years:
        min_year = df_input.index[0]
        max_year = df_input.index[0]
        cut_rows = [*range(min_year, year_start), *range(year_end + 1, max_year + 1)]
        df_input.drop(index=cut_rows, inplace=True)
    return df_input


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
    with_defaults = ["nat_ch4_file", "nat_n2o_file"]
    for key in cfg:
        if key[-4:] == "file":
            if not os.path.exists(cfg[key]):
                if key not in with_defaults:
                    raise FileNotFoundError(f"Input file {key} not found at {cfg[key]}")
                fname = f"natemis_{key.split('_')[1]}.txt"
                cfg[key] = os.path.join(os.getcwd(), "input_OTHER", "NATEMIS", fname)
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
        if "total" not in df_forc.columns:
            if "FORC_NH" in df_forc.columns and "FORC_SH" in df_forc.columns:
                df_forc["total"] = df_forc[["FORC_NH", "FORC_SH"]].mean()
            else:
                df_forc["total"] = df_forc[list(df_forc)].sum(axis=1)
                df_forc["FORC_NH"] = df_forc["total"]
                df_forc["FORC_SH"] = df_forc["total"]
    return df_forc


def read_csv_no_index_col(filename):
    """
    Read input with pandas with no index column

    A method to wrap pandas.read_csv to have index_col = None
    so it can be called as a standard read in method for
    InputHandler.get_data

    Parameters
    ----------
    filename : str
              Path to file to be read in

    Returns
    -------
    pandas.Dataframe
        Dataframe with the intput from the file, possibly cut
        to relevant years
    """
    return pd.read_csv(filename, index_col=None)


class InputHandler:
    """
    Class to handle input sources of different kinds
    and taking care of how to get the correct input
    data
    """

    def __init__(self, cfg):
        """
        Initialise input handler

        Parameters
        ----------
        cfg : dict
           configurations such as input file locations
           Whether run is concentration run, if perturbations
           should be included etc.

        Raises
        ------
        FileNotFoundError
            If forcing file is not found when forcing run is chosen
        """
        self.read_methods = {
            "forc": read_forc,
            "gaspam": read_components,
            "nat_ch4": [read_natural_emissions, "CH4"],
            "nat_n2o": [read_natural_emissions, "N2O"],
            "concentrations": [
                read_inputfile,
                True,
                1750,
                2100,
            ],
            "emissions": self.read_emissions,
            "perturb_em": read_csv_no_index_col,
            "perturb_forc": pd.read_csv,
            "rf_luc": self.read_data_on_year_row,
            "rf_sun": self.read_data_on_year_row,
            "rf_volc_n": self.read_data_on_year_row,
            "rf_volc_s": self.read_data_on_year_row,
        }
        used = {
            "sunvolc": 0,
            "conc_run": False,
            "idtm": 24,
        }
        for key in self.read_methods:
            used[f"{key}_file"] = 0
            used[f"{key}_data"] = 0
        self.cfg = cut_and_check_pamset(
            {"nystart": 1750, "nyend": 2100, "emstart": 1850},
            cfg,
            used=used,
            cut_warnings=True,
        )
        self.read_methods["concentrations"] = [
            read_inputfile,
            True,
            self.cfg["nystart"],
            self.cfg["nyend"],
        ]
        self.set_sun_volc_luc_defaults()

        for nat_emis in ["nat_ch4", "nat_n2o"]:
            if not self.optional_pam(nat_emis) and self.optional_pam("gaspam"):
                self.cfg[f"{nat_emis}_data"] = self.get_flat_natural_emissions(
                    nat_emis.split("_")[1].upper()
                )

        check_inputfiles(self.cfg)

    def get_flat_natural_emissions(self, component):
        """
        Read flat natural emissions from gaspam_file
        if natural emissions are not provided as file or data

        Parameters
        ----------
        component : str
                    String to denote the component, should be
                    N2O or CH4

        Returns
        -------
        pd.Dataframe
                Dataframe of natural emissions for component with years as index
        """
        df_gas = self.get_data("gaspam")
        value = df_gas.at[component, "NAT_EM"]
        years = np.arange(self.cfg["nystart"], self.cfg["nyend"] + 1)
        return pd.DataFrame(data={component: np.ones(len(years)) * value}, index=years)

    def set_sun_volc_luc_defaults(self):
        """
        Set default values for landuse change, solar and
        volcanic data
        """
        if "rf_luc_file" not in self.cfg and "rf_luc_data" not in self.cfg:
            self.cfg["rf_luc_file"] = os.path.join(
                default_data_dir, "IPCC_LUCalbedo.txt"
            )

        if ("sunvolc", 1) in self.cfg.items():
            # Possibly change to allow for other files
            # And for SH to differ from NH
            if "rf_sun_file" not in self.cfg and "rf_sun_data" not in self.cfg:
                self.cfg["rf_sun_file"] = os.path.join(
                    default_data_dir, "solar_IPCC.txt"
                )
            volc_options = [
                "rf_volc_file",
                "rf_volc_data",
                "rf_volc_n_file",
                "rf_volc_s_file",
                "rf_volc_n_data",
                "rf_volc_s_data",
            ]
            intersection = [i for i in volc_options if i in self.cfg]
            if not intersection:
                self.cfg["rf_volc_n_file"] = os.path.join(
                    default_data_dir, "meanVOLCmnd_ipcc_NH.txt"
                )
                self.cfg["rf_volc_s_file"] = os.path.join(
                    default_data_dir, "meanVOLCmnd_ipcc_NH.txt"
                )
            elif "rf_volc_file" in self.cfg:
                self.cfg["rf_volc_n_file"] = self.cfg["rf_volc_file"]
                self.cfg["rf_volc_s_file"] = self.cfg["rf_volc_file"]
            elif "rf_volc_data" in self.cfg:
                self.cfg["rf_volc_n_data"] = self.cfg["rf_volc_data"]
                self.cfg["rf_volc_s_data"] = self.cfg["rf_volc_data"]
        else:
            indices = np.arange(self.cfg["nystart"], self.cfg["nyend"] + 1)
            self.cfg["rf_volc_n_data"] = pd.DataFrame(
                data=np.zeros((self.cfg["nyend"] - self.cfg["nystart"] + 1, 12)),
                index=indices,
                columns=range(12),
            )
            self.cfg["rf_volc_s_data"] = self.cfg["rf_volc_n_data"]
            self.cfg["rf_sun_data"] = pd.DataFrame(
                data=np.zeros(self.cfg["nyend"] - self.cfg["nystart"] + 1),
                index=indices,
            )

    def get_rf_type(self):
        """
        Check if this is rf_run

        Checking if an input key with rf_run info
        is contained in the configuration set, i.e.
        the dataset is an rf_run

        Returns
        -------
        bool
            Whether this is rf_run or not
        """
        rf_type = False
        if "forc_file" in self.cfg or "forc_data" in self.cfg:
            rf_type = True
        return rf_type

    def get_data(self, what):
        """
        Get data from file or sent directly as data

        Method that according to what is looked for
        searches the configuration set self.cfg for
        a parameter that has the same name with _file
        at the end, and if so, looks up the correct
        read in method for that file, or if for
        a parameter with the name with _data at the end,
        in which case it returns that dataset.
        If there is no reading method when name_file
        is present, or if there is no name_file or
        name_data in configurations, then a KeyError
        is raised.

        Parameters
        ----------
        what : str
             what is being looked for, such as emissions,
             concentrations, gaspam etc.

        Returns
        -------
               Dataframe or numpyarray

        Raises
        ------
        KeyError
               If data for what is named is not in configurations
               or if there is no read_method for a file to be read
        """
        if f"{what}_file" in self.cfg:
            if what not in self.read_methods:
                raise KeyError(f"No reading method available for {what}_file")
            if isinstance(self.read_methods[f"{what}"], list):
                return self.read_methods[f"{what}"][0](
                    self.cfg[f"{what}_file"], *self.read_methods[f"{what}"][1:]
                )
            return self.read_methods[f"{what}"](self.cfg[f"{what}_file"])
        if f"{what}_data" in self.cfg:
            # Sanity check of input data?
            return self.cfg[f"{what}_data"]
        raise KeyError(f"No user or default data for {what}")

    def read_emissions(self, filename):
        """
        Read in emission from file

        Method to read in emissions data file by
        calling read_inputfile method and renaming
        CO2 columns

        Parameters
        ----------
        filename : str
             path to emission file to be read in

        Returns
        -------
               pandas.Dataframe
        """
        emis = read_inputfile(
            filename,
            cut_years=True,
            year_start=self.cfg["nystart"],
            year_end=self.cfg["nyend"],
        )
        emis.rename(columns={"CO2": "CO2_FF", "CO2.1": "CO2_AFOLU"}, inplace=True)
        return emis

    def conc_run(self):
        """
        Check if configurations includes conc_run option

        Checking if configurations includes conc_run option
        If included, its value is returned, otherwise False
        is returned

        Returns
        -------
        bool
            Whether this is rf_run or not
        """
        if "conc_run" in self.cfg:
            return self.cfg["conc_run"]
        return False

    def optional_pam(self, which):
        """
        Check if data that can be sent optionally is included

        Checking if an input key for an optionally input data
        is found in the configurations either as file or data
        This could be perturb_em, pertub_forc or forc

        Parameters
        ----------
        which : str
               Type of data to be checked for. This could be
               perturb_em, perturb_forc or forc

        Returns
        -------
        bool
            Whether this is rf_run or not
        """
        return bool(f"{which}_file" in self.cfg or f"{which}_data" in self.cfg)

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

        df_data.set_axis(labels=indices)
        return df_data
