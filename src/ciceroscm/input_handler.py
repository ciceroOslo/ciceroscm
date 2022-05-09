"""
Module to take care of reading in inputs
wither from files or as direct objects
"""
import logging
import os

import numpy as np
import pandas as pd

from ._utils import check_numeric_pamset

LOGGER = logging.getLogger(__name__)


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

                LOGGER.warning(
                    f"Did not find prescribed {key}. Looking in standard path",
                )
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
        self.cfg = check_numeric_pamset(
            {"nystart": 1750, "nyend": 2100, "emstart": 1850}, cfg
        )
        self.read_methods = {
            "forc": read_forc,
            "gaspam": read_components,
            "nat_ch4": [read_natural_emissions, "CH4"],
            "nat_n2o": [read_natural_emissions, "N2O"],
            "concentrations": [
                read_inputfile,
                True,
                self.cfg["nystart"],
                self.cfg["nyend"],
            ],
            "emissions": self.read_emissions,
            "perturb_em": read_csv_no_index_col,
            "perturb_forc": pd.read_csv,
        }
        check_inputfiles(self.cfg)

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
