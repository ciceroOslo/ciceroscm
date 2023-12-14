"""
Perturbation related methods
"""

import numpy as np
import pandas as pd


def calculate_hemispheric_forcing(tracer, q, forc_nh, forc_sh):
    """
    Calculate hemispheric forcing per tracer

    Particular hemisperic factors ar applied for
    SO2, SO4_IND and TROP_O3
    For other compounds the split is equal

    Parameters
    ----------
    tracer : str
          Name of tracer for which to do the split
    q : float
      Global tracer forcing
    forc_nh : float
           Northern hemispheric forcing
    forc_sh : float
           Southern hemispheric forcing

    Returns
    -------
    list
        Containing the two updated hemispheric forcings
        forc_nh, forc_sh
    """
    if tracer in (
        "SO2",
        "SO4_IND",
        "BC",
        "OC",
        "BMB_AEROS",
        "BMB_AEROS_OC",
        "BMB_AEROS_BC",
    ):
        # Smith et al 2020
        forc_nh = forc_nh + q * 1.47
        forc_sh = forc_sh + q * (2 - 1.47)
    elif tracer == "TROP_O3":
        # Skeie et al 2020
        forc_nh = forc_nh + q * 1.45  # 1.3
        forc_sh = forc_sh + q * (2 - 1.45)  # 0.7
    elif tracer == "LANDUSE":
        # Smith et al 2020
        forc_nh = forc_nh + q * 1.42
        forc_sh = forc_sh + q * (2 - 1.42)
    else:
        forc_nh = forc_nh + q
        forc_sh = forc_sh + q
    return forc_nh, forc_sh


def perturb_emissions(input_handler, emissions_df):
    """
    Add emission perturbations to emissions_df

    Read in emissions perturbations file and add
    perturbations to existing emissions

    Parameters
    ----------
    input_handler : ciceroscm.InputHandler
                    InputHandler that takes care of
                    input data
    emissions_df : pandas.Dataframe
                Dataframe of emissions

    Returns
    -------
    pandas.Dataframe
                    Dataframe of emissions with perturbed
                    emissions added in
    """
    pert_df = input_handler.get_data("perturb_em")
    for row in pert_df.itertuples(index=True, name="Pandas"):
        tracer = row.component
        if row.component == "CO2" and "CO2" not in emissions_df:
            tracer = "CO2_FF"
        emissions_df[tracer][row.year] = emissions_df[tracer][row.year] + row.emission
    # might not need to return, as the pandas should change itself
    return emissions_df


class ForcingPerturbation:
    """
    Class to handle  forcing perturbations

    Attributes
    ----------
    perturb_raw : pandas.Dataframe
                  Dataframe storing the perturbations
    years : np.ndarray
            Array of years for which perturbations exist
    compounds : np.ndarray
                Array of components for which perturbations
                exist
    year0 : int
            First year of perturbations
    """

    def __init__(self, input_handler, year0):
        """
        Initialse forcing perturbation instance

        Reading in the perturbations and the years for
        which they apply, also which compounds they apply
        for

        Parameters
        ----------
        perturbation_file : str
                         Path of perturbation file
        year0 : int
             First year of perturbations
        """
        self.perturb_raw = input_handler.get_data("perturb_forc")
        self.years = np.unique(self.perturb_raw["year"].values)
        self.compounds = pd.unique(self.perturb_raw["component"].values)
        self.year0 = year0

    def check_if_year_in_pert(self, year):
        """
        Check if year has perturbations

        Method to check if there are perturbations for a
        given year

        Parameters
        ----------
        year : int
            year to check

        Returns
        -------
        bool
            Truth value of if the year has perturbations

        """
        return bool(year in self.years)

    def check_if_compound_in_pert(self, compound):
        """
        Check if compound is perturbed at some point

        Method to check if a certain compound has perturbations

        Parameters
        ----------
        compound : str
                Name of compound to check

        Returns
        -------
        bool
            Truth value of if the compound has perturbations
        """
        return bool(compound in self.compounds)

    def add_forcing_pert(
        self, totforc, forc_nh, forc_sh, forc, yr
    ):  # pylint: disable=too-many-arguments
        """
        Add forcing perturbations to precalculated forcing

        Method to add perturbations to forcing

        Parameters
        ----------
        totforc : float
               Total forcing for current year
        forc_nh : float
                Northern hemispheric forcing
        forc_sh : float
                Southern hemispheric forcing
        forc : dict
            Of forcing per compound and year
        yr : int
          Year for which to calculate

        Returns
        -------
        list
            Containing updated total and forcing dict and
            hemispheric forcings
            Like this: totforc, forc_nh, forc_sh, forc
        """
        perturb_here = self.perturb_raw[self.perturb_raw.year == yr]
        for row in perturb_here.itertuples(index=True, name="Pandas"):
            value = row.forcing
            tracer = row.component
            totforc = totforc + value
            forc_nh, forc_sh = calculate_hemispheric_forcing(
                tracer, value, forc_nh, forc_sh
            )
            forc[tracer][yr - self.year0] = forc[tracer][yr - self.year0] + value
        return totforc, forc_nh, forc_sh, forc
