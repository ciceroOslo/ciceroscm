"""
Module to handle emission, concentrations and converstion to forcing
"""

import logging
import os

import numpy as np
import pandas as pd

# from ._utils import check_numeric_pamset
from ._utils import cut_and_check_pamset
from .make_plots import plot_output2
from .perturbations import (
    ForcingPerturbation,
    calculate_hemispheric_forcing,
    perturb_emissions,
)
from .pub_utils import make_cl_and_br_dictionaries

LOGGER = logging.getLogger(__name__)


def _rs_function(it, idtm=24):
    """
    Calculate pulse response function for mixed layer

    Calculate pulse response function for mixed layer
    time is the year index*idtm + i, i.e. the month number

    Parameters
    ----------
    it : int
      is the time index, there are idtm time points per yer
    idtm : int
        Number of time points per year, default is 24

    Returns
    -------
    float
         The pulse_response function for this time
    """
    time = it / idtm
    if time < 2.0:
        pulse_response = (
            0.12935
            + 0.21898 * np.exp(-time / 0.034569)
            + 0.17003 * np.exp(-time / 0.26936)
            + 0.24071 * np.exp(-time / 0.96083)
            + 0.24093 * np.exp(-time / 4.9792)
        )
    else:
        pulse_response = (
            0.022936
            + 0.24278 * np.exp(-time / 1.2679)
            + 0.13963 * np.exp(-time / 5.2528)
            + 0.089318 * np.exp(-time / 18.601)
            + 0.03782 * np.exp(-time / 68.736)
            + 0.035549 * np.exp(-time / 232.3)
        )
    return pulse_response


def _rb_function(it, idtm=24):
    """
    Calculate biotic decay function

    Calculate biotic decay function
    time is the year index*idtm + i, i.e. the month number

    Parameters
    ----------
    it : int
      is the time index, there are idtm time points per yer
    idtm : int
        Number of time points per year, default is 24

    Returns
    -------
    float
        The biotic decay function value for this time
    """
    time = it / idtm
    biotic_decay = (
        0.70211 * np.exp(-0.35 * time)
        + 13.4141e-3 * np.exp(-time / 20.0)
        - 0.71846 * np.exp(-55 * time / 120.0)
        + 2.9323e-3 * np.exp(-time / 100.0)
    )
    return biotic_decay


def check_pamset(pamset):
    """
    Check that parameterset has necessary values for run

    Check that parameterset has necessary values for run
    Otherwise set to default values which are defined here

    Parameters
    ----------
    pamset : dict
          Dictionary of parameters to define the physics
          of the run. Values that begin with q are concetration
          or emissions to forcing factors, beta_f is the
          carbon cycle fertilisation factor, and ref_yr is
          the reference year for calculations

    Returns
    -------
    dict
        Updated pamset with default values used where necessary
    """
    required = {
        "qbmb": 0.0,
        "qo3": 0.5,
        "qdirso2": -0.36,
        "qindso2": -0.97,
        "qbc": 0.16,
        "qoc": -0.08,
        "qh2o_ch4": 0.091915,
        "ref_yr": 2010,
        "beta_f": 0.287,
    }

    # pamset = check_numeric_pamset(required, pamset, )
    if "lifetime_mode" not in pamset:
        pamset["lifetime_mode"] = "TAR"

    used = {
        "lifetime_mode": "TAR",
        "just_one": "CO2",
        "idtm": 24,
        "nystart": 1750,
        "nyend": 2100,
        "emstart": 1850,
    }
    pamset = cut_and_check_pamset(required, pamset, used=used, cut_warnings=True)
    return pamset


def check_pamset_consistency(pamset_old, pamset_new):
    """
    Make sure new pamset is consistent with existing instance

    This method is meant to ensure consistency when rerunning
    with the same class instance. Changing to old values
    if they are different for a few values that can't change,
    or setting to old values if none are given

    Parameters
    ----------
    pamset_old : dict
              Original parameterset for instance
    pamset_new : dict
              New pamset for new run
    Returns
    -------
    dict
        The new pamset, augmented to fit the old one
    """
    pams_cant_change = ["idtm", "nystart", "nyend", "emstart", "conc_run"]
    for pam in pams_cant_change:
        if pam in pamset_new and pamset_new[pam] != pamset_old[pam]:
            LOGGER.warning(  # pylint: disable=logging-fstring-interpolation
                f"{pam} can not be changed for same instance of ConcentrationsEmisssionsHandler. Resetting with old value {pamset_old[pam]}. If you want to run with a different value, please create a separate instance",
            )
            pamset_new[pam] = pamset_old[pam]
    for pam, value in pamset_old.items():
        if pam not in pamset_new:
            pamset_new[pam] = value
    return pamset_new


class ConcentrationsEmissionsHandler:
    """
    Class to handle concentrations
    and emissions input for ciceroscm

    Attributes
    ----------
    df_gas : pd.Dataframe
             Dataframe containing names of the various components
             units, possible natural emissions and emissions to
             concentrations and concentrations to forcing conversion
             factors
    conc : dict
           Dictionary with component keys and arrays with concentrations
           for each year as values. These values are calculated
           when emi2conc method is called. Before emissions start
           or if a concentration run is used, they will be read
           directly from conc_in
    nat_emis_ch4 : pd.Dataframe
                   Natural emissions for CH4 per year
    nat_emis_n2o : pd.Dataframe
                   Natural emissions for N2O per year
    pamset : dict
             Dictionary of parameters
    years : np.ndarray
            Array with all the years for the handler
    conc_in : pd.Dataframe
              Input concentrations read from file. These are used
              before emissions start, or throughout the run if
              a concentration run is used
    emis : pd.Dataframe
           Emissions dataframe read from input file
    r_functions : np.ndarray
                  2D array with precalculated values of pulse
                  response and biotic decay functions for each
                  of the idtm values of each year


    """

    # pylint: disable=too-many-instance-attributes
    # Need a few more here or else consider breaking
    # up emissions in it's own class or
    # CO2-handling in its own class

    def __init__(
        self,
        input_handler,
        pamset,
    ):
        """
        Intialising concentrations emissions handler

        Starting by reading in gasparameter file, making
        empty dicts for concentrations and forcing,
        reading in natural emissions, checking the pamset,
        and making an array of years to look at.
        Then reading in concentrations and emissions,
        perturbations if relevant.
        Then initialising pulse and biotic functions,
        finally intialising empty arrays and setting
        CO2 start values.

        Parameters
        ----------
        input_handler : ciceroscm.InputHandler
                       input_handler that takes care of
                       configurations and reading in of
                       data from user
        pamset : dict
           list of physical parameters to define the run
        """
        self.df_gas = input_handler.get_data("gaspam")
        self.conc = {}
        self.forc = {}
        self.nat_emis_ch4 = input_handler.get_data("nat_ch4")
        self.nat_emis_n2o = input_handler.get_data("nat_n2o")
        self.pamset = cut_and_check_pamset(
            {"idtm": 24, "nystart": 1750, "nyend": 2100, "emstart": 1850},
            pamset,
            used={"rs_function": _rs_function, "rb_function": _rb_function},
            cut_warnings=True,
        )
        self.years = np.arange(self.pamset["nystart"], self.pamset["nyend"] + 1)
        self.conc_in = input_handler.get_data("concentrations")
        self.emis = input_handler.get_data("emissions")
        self.pamset["conc_run"] = input_handler.conc_run()
        if input_handler.optional_pam("perturb_em"):
            perturb_emissions(input_handler, self.emis)
        if input_handler.optional_pam("perturb_forc"):
            self.pamset["forc_pert"] = ForcingPerturbation(input_handler, self.years[0])
        self.precalc_r_functions()
        self.pamset["cl_dict"], self.pamset["br_dict"] = make_cl_and_br_dictionaries(
            self.df_gas.index
        )
        # not really needed, but I guess the linter will complain...
        self.reset_with_new_pams(pamset, preexisting=False)

    def reset_with_new_pams(self, pamset, preexisting=True):
        """
        Reset to run again with same emissions etc.

        Resetting arrays and co2 start values also
        making sure parameterset conforms with previously
        existing set if new run is being run

        Parameters
        ----------
        pamset : dict
              Dictionary of physical parameters
        preexisting : bool
                   Defining whether this is the first use of the class
                   or a new run with the same instance, that has to
                   conform to old values
        """
        if preexisting:
            new_pamset = check_pamset(pamset)
            self.pamset = check_pamset_consistency(self.pamset, new_pamset)
        years_tot = len(self.years)
        self.conc = {}
        self.forc = {}
        for tracer in self.df_gas.index:
            if tracer != "CO2.1":
                self.conc[tracer] = {}
                self.forc[tracer] = np.zeros(years_tot)
        self.forc["Total_forcing"] = np.zeros(years_tot)
        self.co2_hold = {
            "yCO2": 0.0,
            "xCO2": 278.0,
            "sCO2": np.zeros(self.pamset["idtm"] * years_tot),
            "emCO2_prev": 0.0,
            "dfnpp": np.zeros(self.pamset["idtm"] * years_tot),
            "ss1": 0.0,
            "sums": 0.0,
        }

    def precalc_r_functions(self):
        """
        Precalculate decay functions either
        sent in pamset or from default

        If functions are sent with keywords rs_function
        or rb_function in the pamset, these must take
        time and number of steps per year as input

        Parameters
        ----------
        pamset
        """
        years_tot = len(self.years)
        self.r_functions = np.empty(
            (2, self.pamset["idtm"] * years_tot)
        )  # if speedup, get this to reflect number of years
        if "rs_function" not in self.pamset:
            self.pamset["rs_function"] = _rs_function
        if "rb_function" not in self.pamset:
            self.pamset["rb_function"] = _rb_function
        self.r_functions[0, :] = [
            self.pamset["rs_function"](it, self.pamset["idtm"])
            for it in range(self.pamset["idtm"] * years_tot)
        ]
        self.r_functions[1, :] = [
            self.pamset["rb_function"](it, self.pamset["idtm"])
            for it in range(self.pamset["idtm"] * years_tot)
        ]

    def calculate_strat_quantities(self, yr):
        """
        Calculate sumcl and sumbr in stratosphere

        Calculating sum of stratospheric effects of
        chlorinated gase and sum of stratospheric effects
        of halon quantities.

        Parameters
        ----------
        yr : int
          Year for which to calculate

        Returns
        -------
        list
            Containing the sumcl, value for chlorinated gases
            and the sumbr value for halon gases.
        """
        yr0 = int(self.years[0])
        if yr <= yr0:
            return 0.0, 0.0
        sumcl = 0
        sumbr = 0
        for comp, mult in self.pamset["cl_dict"].items():
            num = self.conc[comp][yr] - self.conc[comp][yr0]
            # sumcl = sumcl + (mult * (self.conc[comp][yr] - self.conc[comp][yr0])) ** 1.7
            sumcl = sumcl + np.sign(num) * np.power(mult * abs(num), 1.7)
        for comp, mult in self.pamset["br_dict"].items():
            sumbr = sumbr + mult * (self.conc[comp][yr] - self.conc[comp][yr0])
        return sumcl, sumbr

    def calculate_forc_three_main(self, yr):  # pylint: disable=too-many-locals
        """
        Calculate forcings for components CO2, N2O and CH4

        Method to calculate the forcing from the main three componenents
        as these three are a bit intertwined. Using a coefficient matrix

        Parameters
        ----------
        yr : int
          Year for which to calculate

        Returns
        -------
        list
            Containing the total forcing from the three, and the
            combinded forcing value on each hemisphere
            tot_forc, forc_nh, forc_sh
        """
        # Setting some constants. What are they? Should they be parametrisable?
        # (Etminan I guess...)

        # Array with coefficients
        co2_n2o_ch4_coeff = np.array(
            [
                [-2.4e-7, 7.2e-4, -2.1e-4],
                [-8.0e-6, 4.2e-6, -4.9e-6],
                [-1.3e-6, -8.2e-6, None],
            ]
        )
        # a1 = -2.4e-7
        # b1 = 7.2e-4
        # c1 = -2.1e-4
        # a2 = -8.0e-6
        # b2 = 4.2e-6
        # c2 = -4.9e-6
        # a3 = -1.3e-6
        # b3 = -8.2e-6
        yr_0 = self.years[0]
        c0_n2o = self.conc["N2O"][yr_0]
        c_n2o = self.conc["N2O"][yr]
        c0_co2 = self.conc["CO2"][yr_0]
        c_co2 = self.conc["CO2"][yr]
        c0_ch4 = self.conc["CH4"][yr_0]
        c_ch4 = self.conc["CH4"][yr]
        q_co2 = (
            co2_n2o_ch4_coeff[0, 0] * (c_co2 - c0_co2) ** 2
            + co2_n2o_ch4_coeff[0, 1] * np.abs(c_co2 - c0_co2)
            + co2_n2o_ch4_coeff[0, 2] * 0.5 * (c0_n2o + c_n2o)
            + 5.36
        ) * np.log(
            c_co2 / c0_co2
        )  # Etminan et al 2016 RF

        # Note there were a few unused calculations in the original
        # for n2o and ch4 including som fmn values
        q_n2o = (
            co2_n2o_ch4_coeff[1, 0] * 0.5 * (c_co2 + c0_co2)
            + co2_n2o_ch4_coeff[1, 1] * 0.5 * (c_n2o + c0_n2o)
            + co2_n2o_ch4_coeff[1, 2] * 0.5 * (c_ch4 + c0_ch4)
            + 0.117
        ) * (np.sqrt(c_n2o) - np.sqrt(c0_n2o))
        q_ch4 = (
            co2_n2o_ch4_coeff[2, 0] * 0.5 * (c_ch4 + c0_ch4)
            + co2_n2o_ch4_coeff[2, 1] * 0.5 * (c_n2o + c0_n2o)
            + 0.043
        ) * (np.sqrt(c_ch4) - np.sqrt(c0_ch4))
        # Feedback factor: Smith et al 2018
        q_co2 = self.df_gas["SARF_TO_ERF"]["CO2"] * q_co2  # + FORC_PERT(yr_ix,trc_ix))
        q_n2o = self.df_gas["SARF_TO_ERF"]["N2O"] * q_n2o  # + FORC_PERT(yr_ix,trc_ix))
        q_ch4 = self.df_gas["SARF_TO_ERF"]["CH4"] * q_ch4  # + FORC_PERT(yr_ix,trc_ix))

        self.forc["CO2"][yr - yr_0] = q_co2
        self.forc["CH4"][yr - yr_0] = q_ch4
        self.forc["N2O"][yr - yr_0] = q_n2o
        tot_forc = q_co2 + q_n2o + q_ch4
        forc_nh, forc_sh = calculate_hemispheric_forcing(
            "three_main", tot_forc, 0.0, 0.0
        )

        return tot_forc, forc_nh, forc_sh

    def tropospheric_ozone_forcing(self, yr):
        """
        Calculate tropospheric ozone forcing

        Using the emissions of various other gases
        Before emissions start, the fossil fuel CO2
        change from the start is used
        After emissions start, concentrations of
        methane are combined with emissions of NOx, CO and NMVOC
        to calculate a total troposperic ozone forcing

        Parameters
        ----------
        yr : int
          Year for which to calculate

        Returns
        -------
        float
             tropospheric ozone forcing
        """
        emstart = self.pamset["emstart"]
        yr_0 = self.years[0]
        yr_emstart = emstart - yr_0
        yr_ix = yr - yr_0
        tracer = "TROP_O3"
        if yr_ix < yr_emstart:
            # Uses change in CO2_FF emissions
            if self.emis["CO2_FF"][self.pamset["ref_yr"]] != self.emis["CO2_FF"][yr_0]:
                q = (
                    (self.emis["CO2_FF"][yr] - self.emis["CO2_FF"][yr_0])
                    / (
                        self.emis["CO2_FF"][self.pamset["ref_yr"]]
                        - self.emis["CO2_FF"][yr_0]
                    )
                    * self.pamset["qo3"]
                )
            else:
                q = (
                    (self.emis["CO2_FF"][yr] - self.emis["CO2_FF"][yr_0])
                ) * self.pamset["qo3"]

        else:
            # ALOG(1700.0))  !Concentration in 2010 &
            self.conc[tracer][yr] = (
                30.0
                + 6.7
                * (
                    np.log(self.conc["CH4"][yr])
                    - np.log(self.conc_in["CH4"][self.pamset["ref_yr"]])
                )
                + 0.17
                * (self.emis["NOx"][yr] - self.emis["NOx"][self.pamset["ref_yr"]])
                + 0.0014
                * (self.emis["CO"][yr] - self.emis["CO"][self.pamset["ref_yr"]])
                + 0.0042
                * (self.emis["NMVOC"][yr] - self.emis["NMVOC"][self.pamset["ref_yr"]])
            )
            # RBS101115
            # IF (yr_ix.LT.yr_2010) THEN ! Proportional to TROP_O3 build-up
            # Rewritten a bit, place to check for differences...
            forc_pre_emstart = self.forc[tracer][yr_emstart - 1]
            value = self.conc[tracer][yr]
            value_0 = self.conc[tracer][self.pamset["emstart"]]
            q = forc_pre_emstart + (value - value_0) / (30.0 - value_0) * (
                self.pamset["qo3"] - forc_pre_emstart
            )
        return q

    def conc2forc(self, yr, rf_luc, rf_sun):  # pylint: disable=too-many-branches
        """
        Calculate forcing from concentrations

        Looping through all species and getting the
        forcing from their concentrations. First calling
        a method to calculate intial forcing for three main
        species, CO2, N2O and CH4. This will also initialise
        total forcing and hemispherically split forcing.
        Then looping over species
        calculated dirctly from emissions (SO2, SO4_IND, OC and BC)
        Then looping over other gases with alpha factor,
        then ozone in troposphere and stratosphere, and stratospheric water
        vapour, finally other species, and adding perturbations if
        present. For every species add to total forcing and hemispherically
        split forcings, finally add solar forcing.

        Parameters
        ----------
        yr : int
          Year fro which to calculate
        rf_luc : float
              Land use change forcing
        rf_sun : float
              Solar forcing

        Returns
        -------
        list
            Consisting of total forcing and hemispherically split
            forcing for the year.
            In this way: tot_forc, forc_nh, forc_sh
        """
        ref_emission_species = {
            "SO2": ["SO2", self.pamset["qdirso2"]],
            "SO4_IND": ["SO2", self.pamset["qindso2"]],
            "OC": ["OC", self.pamset["qoc"]],
            "BC": ["BC", self.pamset["qbc"]],
            "BMB_AEROS": ["BMB_AEROS_OC", self.pamset["qbmb"]],
        }
        # Intialising with the combined values from CO2, N2O and CH4
        tot_forc, forc_nh, forc_sh = self.calculate_forc_three_main(yr)
        yr_0 = self.years[0]
        # Finish per tracer calculations, add per tracer to printable df and sum the total
        for tracer in self.df_gas.index:
            if tracer in ["CO2", "N2O", "CH4"]:
                continue
            q = 0
            if tracer == "LANDUSE":
                q = rf_luc
            elif tracer in ref_emission_species:
                # Natural emissions
                # (after IPCC TPII on simple climate models, 1997)
                # enat = 42.0 Not used, why is this here?
                # Emission in reference year
                # SO2, SO4_IND, BC and OC are treated exactly the same
                # Only with emission to concentration factors differing
                # These are held in dictionary
                erefyr = (
                    self.emis[ref_emission_species[tracer][0]][self.pamset["ref_yr"]]
                    - self.emis[ref_emission_species[tracer][0]][yr_0]
                )
                if erefyr != 0.0:  # pylint: disable=compare-to-zero
                    frac_em = (
                        self.emis[ref_emission_species[tracer][0]][yr]
                        - self.emis[ref_emission_species[tracer][0]][yr_0]
                    ) / erefyr
                    q = ref_emission_species[tracer][1] * frac_em

            elif (
                tracer in self.df_gas.index
                and self.df_gas["ALPHA"][tracer] != 0  # pylint: disable=compare-to-zero
            ):
                q = (
                    (self.conc[tracer][yr] - self.conc[tracer][yr_0])
                    * self.df_gas["ALPHA"][tracer]
                    * self.df_gas["SARF_TO_ERF"][tracer]
                )  # +forc_pert
            elif tracer == "TROP_O3":
                q = (
                    self.tropospheric_ozone_forcing(yr)
                    * self.df_gas["SARF_TO_ERF"][tracer]
                )
            elif tracer == "STRAT_O3":
                sumcl, sumbr = self.calculate_strat_quantities(yr - 3)
                # Updated according to IPCC AR4.
                # Multiply by factor 0.287737 (=-0.05/-0.17377, AR4/SCM)
                q = (
                    -(
                        # self.pamset["qo3"]
                        0.05
                        / 0.17377
                        * (0.000552 * (sumcl) + 3.048 * sumbr)
                    )
                    / 1000.0
                    * self.df_gas["SARF_TO_ERF"][tracer]
                )
            elif tracer == "STRAT_H2O":
                q = (
                    self.pamset["qh2o_ch4"] * self.forc["CH4"][yr - yr_0]
                ) * self.df_gas["SARF_TO_ERF"][
                    tracer
                ]  # + FORC_PERT(yr_ix,trc_ix)
            elif tracer == "OTHER":
                # Possible with forcing perturbations for other
                # components such as contrails, cirrus etc...
                pass
            self.forc[tracer][yr - yr_0] = q  # + FORC_PERT(yr_ix,trc_ix)
            # Calculating hemispheric forcing:
            forc_nh, forc_sh = calculate_hemispheric_forcing(
                tracer, q, forc_nh, forc_sh
            )
            tot_forc = tot_forc + q
            # print("Forcer: %s, tot_forc: %f, FN: %f, FS: %f, q: %f"%(tracer, tot_forc, forc_nh, forc_sh, q)

        # Adding forcing perturbations if they exist:
        if "forc_pert" in self.pamset:
            if self.pamset["forc_pert"].check_if_year_in_pert(yr):
                tot_forc, forc_nh, forc_sh, self.forc = self.pamset[
                    "forc_pert"
                ].add_forcing_pert(tot_forc, forc_nh, forc_sh, self.forc, yr)

        # Adding solar forcing
        # tot_forc = tot_forc + rf_sun
        if "just_one" in self.pamset:
            tot_forc = self.forc[self.pamset["just_one"]][yr - yr_0]
            forc_nh, forc_sh = calculate_hemispheric_forcing(
                self.pamset["just_one"], tot_forc, 0, 0
            )
        self.forc["Total_forcing"][yr - yr_0] = tot_forc
        forc_nh = forc_nh + rf_sun
        forc_sh = forc_sh + rf_sun

        return tot_forc, forc_nh, forc_sh

    def emi2conc(self, yr):
        """
        Calculate concentrations from emissions

        If conc_run or prior to emissions start, concentrations
        are simply read in from existing data frame, after that
        all tracers are looped over to calculate concentrations
        from emissions. For CO2 this is doen in a separate method.
        For CH4 and N2O, natural emissions are added.

        Parameters
        ----------
        yr : int
          Year for which to calculate
        """
        # Do per tracer emissions to concentrations, update concentrations df
        # NBNB! Remember to move calculation of  Trop_O3 concentration
        # here from conc2forc.
        # First calculate O3 concentrations (AS IN TAR p269, table 4.11 B)
        # CONC(yr_ix,trc_ix) =  30.0 + 6.7 * &
        # (ALOG(CONC(yr_ix,trcID("CH4")))-ALOG(1832.0)) &   !ALOG(1700.0))  !Concentration in 2010 &
        #  + 0.17 * (EMISSIONS(yr_ix,trcID("NOx"))-EM2010(trcID("NOx"))) &
        #  + 0.0014 * (EMISSIONS(yr_ix,trcID("CO"))-EM2010(trcID("CO"))) &
        #  + 0.0042 *(EMISSIONS(yr_ix,trcID("NMVOC"))-EM2010(trcID("NMVOC")))
        if self.pamset["conc_run"]:
            self.fill_one_row_conc(yr)
            return
        # Before emissions start
        if yr < self.pamset["emstart"]:
            self.co2em2conc(yr)
            self.fill_one_row_conc(yr, avoid=["CO2"])
            return
        self.add_row_of_zeros_conc(yr)

        for tracer in self.df_gas.index:
            # something something postscenario...?
            if self.df_gas["CONC_UNIT"][tracer] == "-":
                # Forcing calculated from emissions
                continue
            if tracer == "CO2":
                self.co2em2conc(yr)
                continue
            if yr < self.pamset["emstart"]:
                self.fill_one_row_conc(yr)

            if yr > self.years[0]:
                conc_local = self.conc[tracer][yr - 1]
            else:
                conc_local = self.conc_in[tracer][yr]

            q = 1.0 / self.df_gas["TAU1"][tracer]

            if tracer == "CH4":
                self.df_gas.at[tracer, "NAT_EM"] = self.nat_emis_ch4["CH4"][yr]
                q = self.methane_lifetime(q, conc_local, yr)
            if tracer == "N2O":
                self.df_gas.at[tracer, "NAT_EM"] = self.nat_emis_n2o["N2O"][yr]

            # natural emissions, from gasspamfile
            emis = self.emis[tracer][yr] + self.df_gas["NAT_EM"][tracer]

            point_conc = emis / self.df_gas["BETA"][tracer]
            # Rewrote this quite a bit from an original loop,
            # but I think it is mathematically equivalent
            self.conc[tracer][yr] = point_conc / q + (
                conc_local - point_conc / q
            ) * np.exp(-q)

    def methane_lifetime(self, q, conc_local, yr):
        """
        Calculate methane concentrations from emissions

        Calculate methane concentrations from emissions
        using different liftetime modes to calculate

        Parameters
        ----------
        q : float
         Forcing without feedbacks
        conc_local : float
                  Concentration of methane in previous timestep
        yr : int
          Year for which to calculate
        Returns
        -------
        float
             Concentrations adjusted for lifetime / feedback
        """
        ch4_wigley_exp = -0.238
        if self.pamset["lifetime_mode"] == "TAR":
            # 1751 is reference conc in 2000
            dln_oh = (
                -0.32 * (np.log(conc_local) - np.log(1751.0))
                + 0.0042 * (self.emis["NOx"][yr] - self.emis["NOx"][2000])
                - 0.000105 * (self.emis["CO"][yr] - self.emis["CO"][2000])
                - 0.000315 * (self.emis["NMVOC"][yr] - self.emis["NMVOC"][2000])
            )
            q = q * (dln_oh + 1)
        elif self.pamset["lifetime_mode"] == "CONSTANT_12":
            q = 1.0 / 12.0
        elif self.pamset["lifetime_mode"] == "WIGLEY":
            q = q * (((conc_local / 1700.0)) ** (ch4_wigley_exp))

        q = q + 1.0 / self.df_gas["TAU2"]["CH4"] + 1.0 / self.df_gas["TAU3"]["CH4"]

        return q

    def co2em2conc(self, yr):  # pylint: disable=too-many-locals
        """
        Calculate co2 concentrations from emissions

        Method to calculate co2 concentrations from emissions
        Implementing a rudimentary carbon cycle which loops over
        idtm (usually 24) timesteps a year

        Parameters
        ----------
        yr : int
          Year for which to calculate
        """
        # Area of the ocean (m^2)
        ocean_area = 3.62e14

        # Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
        coeff = 1.0 / (ocean_area * 9.06)

        # TIMESTEP (YR)
        dt = 1.0 / self.pamset["idtm"]

        # Conversion factor ppm/kg --> umol*/m3
        conv_factor = 1.722e17

        # USING MIXED LAYER DEPTH = 75 metres
        mixed_layer_depth = 75.0

        cc1 = dt * ocean_area * coeff / (1 + dt * ocean_area * coeff / 2.0)
        yr_ix = yr - self.years[0]
        # Monthloop:
        em_co2_common = (
            self.emis["CO2_FF"][yr]
            + self.emis["CO2_AFOLU"][yr]
            + self.df_gas["NAT_EM"]["CO2"]
        )
        for i in range(self.pamset["idtm"]):
            it = yr_ix * self.pamset["idtm"] + i
            sumf = 0.0

            # Net emissions, including biogenic fertilization effects
            if it > 0:
                self.co2_hold["dfnpp"][it] = (
                    60 * self.pamset["beta_f"] * np.log(self.co2_hold["xCO2"] / 278.0)
                )
            if it > 0:
                sumf = float(
                    np.dot(
                        self.co2_hold["dfnpp"][1:it],
                        np.flip(self.r_functions[1, : it - 1]),
                    )
                )

            ffer = self.co2_hold["dfnpp"][it] - dt * sumf
            em_co2 = (em_co2_common - ffer) / 2.123

            if it == 0:  # pylint: disable=compare-to-zero
                self.co2_hold["ss1"] = 0.5 * em_co2 / (ocean_area * coeff)
                ss2 = self.co2_hold["ss1"]
                self.co2_hold["sums"] = 0.0
            else:
                ss2 = 0.5 * em_co2 / (ocean_area * coeff) - self.co2_hold["yCO2"] / (
                    dt * ocean_area * coeff
                )
                self.co2_hold["sums"] = (
                    self.co2_hold["sums"]
                    + self.co2_hold["emCO2_prev"] / (ocean_area * coeff)
                    - self.co2_hold["sCO2"][it - 1]
                )
            self.co2_hold["sCO2"][it] = cc1 * (
                self.co2_hold["sums"] + self.co2_hold["ss1"] + ss2
            )
            self.co2_hold["emCO2_prev"] = em_co2
            if it > 0:
                sumz = np.dot(
                    self.co2_hold["sCO2"][: it - 1], np.flip(self.r_functions[0, 1:it])
                )
            else:
                sumz = 0.0

            z_co2 = (
                conv_factor
                * coeff
                * dt
                / mixed_layer_depth
                * (sumz + 0.5 * self.co2_hold["sCO2"][it])
            )
            self.co2_hold["yCO2"] = (
                1.3021 * z_co2
                + 3.7929e-3 * (z_co2**2)
                + 9.1193e-6 * (z_co2**3)
                + 1.488e-8 * (z_co2**4)
                + 1.2425e-10 * (z_co2**5)
            )
            self.co2_hold["xCO2"] = (
                self.co2_hold["sCO2"][it] + self.co2_hold["yCO2"] + 278.0
            )
            # print("it: %d, emCO2: %e, sCO2: %e, zCO2: %e, yCO2: %e, xCO2: %e, ss1: %e, ss2: %e, dnfpp:%e"%(it, em_co2, self.co2_hold["sCO2"][it], z_co2, self.co2_hold["yCO2"], self.co2_hold["xCO2"], self.co2_hold["ss1"], ss2, self.co2_hold["dfnpp"][it]))
        self.conc["CO2"][yr] = self.co2_hold["xCO2"]

    def fill_one_row_conc(self, yr, avoid=None):
        """
        Fill in one row of concentrations in conc_dict

        Fill in a row of zeros in concentrations dictionary
        from prescribed concentrations. If some traces are
        to be avoided (typically CO2 calculated from emissions)
        a list of them can be sent as inputs

        Parameters
        ----------
        yr : int
          Year for which to add concentrations from prescribed
        avoid : list
             optional list of tracers for which not to read
             presecribed concentrations. Typically CO2 when CO2
             is calculated from emissions, while other compounds
             are read from prescribed data.
        """
        for tracer, value_dict in self.conc.items():
            if avoid and tracer in avoid:
                continue
            if tracer in self.conc_in:
                value_dict[yr] = self.conc_in[tracer][yr]
            else:
                value_dict[yr] = 0

    def add_row_of_zeros_conc(self, yr):
        """
        Fill in one row of concentrations in conc_dict

        Fill in a row of zeros in concentrations dictionary
        for a single year

        Parameters
        ----------
        yr : int
          Year for which to add zeros
        """
        for value_dict in self.conc.values():
            value_dict[yr] = 0

    def write_output_to_files(self, cfg, make_plot=False):
        """
        Write results to files after run

        Write results for emissions, concentrations
        and forcings to files. Format is as in the
        fortran implementation, with separate files for
        emissions (emis), concentrations (conc) and
        forcing (forc)

        Parameters
        ----------
        cfg : dict
           Configurations to define where to put output
           files and what prefix to have for file name
        make_plot : bool
           Whether the output should be plottet or not
        """
        if "output_folder" in cfg:
            # Make os independent?
            outdir = os.path.join(os.getcwd(), cfg["output_folder"])
        else:
            outdir = os.path.join(os.getcwd(), "output")

        df_forc = pd.DataFrame(data=self.forc, index=self.years)
        df_forc["Year"] = self.years

        cols = df_forc.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_forc = df_forc[cols]
        df_forc.rename(columns={"SO2": "SO4_DIR"}, inplace=True)
        df_emis = self.emis.drop(labels=["CO2_FF", "CO2_AFOLU"], axis=1).drop(
            labels=np.arange(self.years[-1] + 1, self.emis.index[-1] + 1), axis=0
        )
        df_emis["Year"] = self.years
        df_emis["CO2"] = self.emis["CO2_FF"] + self.emis["CO2_AFOLU"]
        cols = df_emis.columns.tolist()
        cols = cols[-2:] + cols[:-2]
        df_emis = df_emis[cols]
        self.conc["Year"] = self.years

        df_conc = pd.DataFrame(data=self.conc, index=self.years)
        cols = df_conc.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_conc = df_conc[cols]
        for tracer in self.df_gas.index:
            if tracer not in df_emis.columns.tolist():
                df_emis[tracer] = np.zeros(len(self.years))
        if "output_prefix" in cfg:
            filename_start = cfg["output_prefix"]
        else:
            filename_start = "output"
        df_forc.to_csv(
            os.path.join(outdir, f"{filename_start}_forc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_conc.to_csv(
            os.path.join(outdir, f"{filename_start}_conc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )

        df_emis.to_csv(
            os.path.join(outdir, f"{filename_start}_em.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )

        if make_plot:
            plot_output2("forc", df_forc, outdir)
            plot_output2("emis", df_emis, outdir, self.df_gas["EM_UNIT"])
            plot_output2("conc", df_conc, outdir, self.df_gas["CONC_UNIT"])

    def add_results_to_dict(self):
        """
        Adding results to results dictionary

        Returns
        -------
        dict
            Containing run emissions, concentrations and forcings
            in the form of pandas.Dataframe with years and tracers
        """
        df_forc = pd.DataFrame(data=self.forc, index=self.years)
        df_forc["Year"] = self.years

        cols = df_forc.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_forc = df_forc[cols]
        df_forc.rename(columns={"SO2": "SO4_DIR"}, inplace=True)
        df_emis = self.emis.drop(labels=["CO2_FF", "CO2_AFOLU"], axis=1).drop(
            labels=np.arange(self.years[-1] + 1, self.emis.index[-1] + 1), axis=0
        )
        df_emis["Year"] = self.years
        df_emis["CO2"] = self.emis["CO2_FF"] + self.emis["CO2_AFOLU"]
        cols = df_emis.columns.tolist()
        cols = cols[-2:] + cols[:-2]
        df_emis = df_emis[cols]
        self.conc["Year"] = self.years

        df_conc = pd.DataFrame(data=self.conc, index=self.years)
        cols = df_conc.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_conc = df_conc[cols]
        results = {}
        results["emissions"] = df_emis
        results["concentrations"] = df_conc
        results["forcing"] = df_forc
        return results
