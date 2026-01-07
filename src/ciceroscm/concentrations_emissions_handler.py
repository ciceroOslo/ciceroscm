"""
Module to handle emission, concentrations and converstion to forcing
"""

import logging
import os

import numpy as np
import pandas as pd

# from ._utils import check_numeric_pamset
from ._utils import cut_and_check_pamset
from .carbon_cycle.common_carbon_cycle_functions import calculate_airborne_fraction
from .component_factory_functions import create_carbon_cycle_model
from .make_plots import plot_output2
from .perturbations import (
    ForcingPerturbation,
    calculate_hemispheric_forcing,
    perturb_emissions,
)
from .pub_utils import make_cl_and_br_dictionaries

LOGGER = logging.getLogger(__name__)


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
          or emissions to forcing factors
          and ref_yr is the reference year for calculations

    Returns
    -------
    dict
        Updated pamset with default values used where necessary
    """
    if pamset is None:
        pamset = {}
    required = {
        "qbmb": 0.0,
        "qo3": 0.5,
        "qdirso2": -0.00308,
        "qindso2": -0.97 / 57.052577209999995,
        "qbc": 0.0279,
        "qoc": -0.00433,
        "qh2o_ch4": 0.091915,
        "ref_yr": 2010,
        "qnmvoc": 0.0,
        "qnh3": 0.0,
        "qnox": 0.0,
        "idtm": 24,
        "nystart": 1750,
        "nyend": 2100,
        "emstart": 1850,
    }

    # pamset = check_numeric_pamset(required, pamset, )
    if "lifetime_mode" not in pamset:
        pamset["lifetime_mode"] = "TAR"
    if "carbon_cycle_model" not in pamset:
        pamset["carbon_cycle_model"] = "default"
    used = {
        "lifetime_mode": "TAR",
        "just_one": "CO2",
        "carbon_cycle_model": "default",
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

    def __init__(self, input_handler, pamset, pamset_carbon=None):
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
        self.pamset = check_pamset(pamset)
        self.years = np.arange(self.pamset["nystart"], self.pamset["nyend"] + 1)
        self.conc_in = input_handler.get_data("concentrations")
        self.emis = input_handler.get_data("emissions")
        self.pamset["conc_run"] = input_handler.conc_run()
        if input_handler.optional_pam("perturb_em"):
            perturb_emissions(input_handler, self.emis)
        if input_handler.optional_pam("perturb_forc"):
            self.pamset["forc_pert"] = ForcingPerturbation(input_handler, self.years[0])
        self.pamset["cl_dict"], self.pamset["br_dict"] = make_cl_and_br_dictionaries(
            self.df_gas.index
        )
        # Precalculating for vanilla gases
        self.precalc_dict = {}
        self._precalculate_vanilla_gases()

        # Setting up carbon cycle model
        model_type = self.pamset["carbon_cycle_model"]  # Default to "default"
        self.carbon_cycle = create_carbon_cycle_model(
            model_type, self.pamset, pamset_carbon
        )
        # not really needed, but I guess the linter will complain...
        self.reset_with_new_pams(pamset, pamset_carbon, preexisting=False)

    def _precalculate_vanilla_gases(self):
        df_gases_vanilla = self.df_gas.copy()

        to_drop = ["CO2", "CH4", "N2O"]
        for index, row in df_gases_vanilla.iterrows():
            if row["CONC_UNIT"] == "-":
                to_drop.append(index)
        df_gases_vanilla.drop(labels=to_drop, inplace=True)
        self.precalc_dict["df_gas_vanilla"] = df_gases_vanilla
        self.precalc_dict["precalc_conc"] = (
            self._precalculate_concentrations_vanilla_gases(df_gases_vanilla)
        )
        self.precalc_dict["precalc_erf"] = self._precalculate_erf_vanilla_gases()

        self._add_precalculated_strat_o3()

        if "forc_pert" in self.pamset and self.pamset["forc_pert"] is not None:
            self.pamset["forc_pert"].add_forcing_pert_for_vanilla(
                self.precalc_dict["precalc_erf"]
            )

        self.precalc_dict["precalc_erf"]["TOT_vanilla"] = self.precalc_dict[
            "precalc_erf"
        ].sum(axis=1)
        self.precalc_dict["gases_list_spicy"] = []
        for index in self.df_gas.index:
            if index not in self.precalc_dict["precalc_erf"].columns:
                self.precalc_dict["gases_list_spicy"].append(index)

    def avoid_regional_double_counting(self):
        """
        Set global aerosol focing parameters to zero if regional forcing is used to avoid double counting
        """
        aerosols = {
            "SO2": "qdirso2",
            "BC": "qbc",
            "OC": "qoc",
        }
        for aerosol, aerforc in aerosols.items():
            if any(comp.startswith(f"{aerosol}_") for comp in self.df_gas.index):
                self.pamset[aerforc] = 0.0

    def _precalculate_erf_vanilla_gases(self):
        """
        Precalculate erf for vanilla gases
        """
        conc0 = self.precalc_dict["precalc_conc"].iloc[0, :].to_numpy()
        alpha = self.precalc_dict["df_gas_vanilla"]["ALPHA"].to_numpy()
        sarf_to_erf = self.precalc_dict["df_gas_vanilla"]["SARF_TO_ERF"].to_numpy()

        erf_data = (
            sarf_to_erf * alpha * (self.precalc_dict["precalc_conc"].to_numpy() - conc0)
        )

        erf_df = pd.DataFrame(
            data=erf_data,
            columns=self.precalc_dict["precalc_conc"].columns,
            index=self.precalc_dict["precalc_conc"].index,
        )
        return erf_df

    def _precalculate_concentrations_vanilla_gases(self, df_gases_vanilla):
        conc_in_vanilla = self.conc_in.copy()[df_gases_vanilla.index]
        if self.pamset["conc_run"] or self.pamset["nyend"] <= self.pamset["emstart"]:
            return conc_in_vanilla.iloc[
                self.pamset["nystart"]
                - conc_in_vanilla.index[0] : self.pamset["nyend"]
                - conc_in_vanilla.index[0]
                + 1
            ]

        conc_in_vanilla = conc_in_vanilla.iloc[
            self.pamset["nystart"]
            - conc_in_vanilla.index[0] : self.pamset["emstart"]
            - conc_in_vanilla.index[0]
        ]

        emis_vanilla = self.emis.copy()[df_gases_vanilla.index]
        emis_vanilla = emis_vanilla.iloc[
            self.pamset["emstart"]
            - emis_vanilla.index[0] : self.pamset["nyend"]
            - emis_vanilla.index[0]
            + 1
        ]
        q = 1.0 / df_gases_vanilla["TAU1"].to_numpy()
        emis = emis_vanilla.to_numpy() + df_gases_vanilla["NAT_EM"].to_numpy()
        conc_rows = []
        for i in range(emis.shape[0]):
            emis_now = emis[i, :]
            if i == 0:
                conc_prev = conc_in_vanilla.iloc[-1, :].to_numpy()
            else:
                conc_prev = conc_rows[-1]
            conc = emis_now / df_gases_vanilla["BETA"].to_numpy() / q * (
                1 - np.exp(-q)
            ) + conc_prev * np.exp(-q)
            conc_rows.append(conc)
        conc_rows_df = pd.DataFrame(
            data=np.array(conc_rows),
            columns=conc_in_vanilla.columns,
            index=emis_vanilla.index,
        )
        return pd.concat([conc_in_vanilla, conc_rows_df])

    def _add_precalculated_strat_o3(self):
        sumcl, sumbr = np.zeros(len(self.years)), np.zeros(len(self.years))
        sumcl[3:], sumbr[3:] = self.calculate_strat_quantities(
            self.years[:-3], self.precalc_dict["precalc_conc"]
        )
        q = (
            -(
                # self.pamset["qo3"]
                0.05
                / 0.17377
                * (0.000552 * (sumcl) + 3.048 * sumbr)
            )
            / 1000.0
            * self.df_gas["SARF_TO_ERF"]["STRAT_O3"]
        )
        self.precalc_dict["precalc_erf"]["STRAT_O3"] = q

    def reset_with_new_pams(self, pamset, pamset_carbon=None, preexisting=True):
        """
        Reset to run again with same emissions etc.

        Resetting arrays and co2 start values also
        making sure parameterset conforms with previously
        existing set if new run is being run

        Parameters
        ----------
        pamset : dict
              Dictionary of physical parameters
        pamset_carbon : dict
            Dictionary of parameters to be passed to the carbon cycle model
        preexisting : bool
                   Defining whether this is the first use of the class
                   or a new run with the same instance, that has to
                   conform to old values
        """
        if preexisting:
            new_pamset = check_pamset(pamset)
            self.pamset = check_pamset_consistency(self.pamset, new_pamset)

        self.avoid_regional_double_counting()
        self.carbon_cycle.reset_co2_hold(pamset_carbon)

        years_tot = len(self.years)
        self.conc = {}
        self.forc = {}
        for tracer in self.precalc_dict["gases_list_spicy"]:
            if tracer != "CO2.1":
                if tracer == "TROP_O3" or self.df_gas["CONC_UNIT"][tracer] != "-":
                    self.conc[tracer] = {}
                self.forc[tracer] = np.zeros(years_tot)
        self.forc["Total_forcing"] = np.zeros(years_tot)

    def calculate_strat_quantities(self, yr, conc):
        """
        Calculate sumcl and sumbr in stratosphere

        Calculating sum of stratospheric effects of
        chlorinated gase and sum of stratospheric effects
        of halon quantities.

        Parameters
        ----------
        yr : int or np.ndarray
          Year for which to calculate
        conc : dict or pd.DataFrame
            Concetrations DataFrame or dict from which
            concentrations of chlorinated and halon gases
            for the year or years can be read

        Returns
        -------
        list
            Containing the sumcl, value for chlorinated gases
            and the sumbr value for halon gases.
        """
        yr0 = int(self.years[0])
        if np.isscalar(yr):
            if yr <= yr0:
                return 0.0, 0.0
        sumcl = 0
        sumbr = 0
        for comp, mult in self.pamset["cl_dict"].items():
            num = conc[comp][yr] - conc[comp][yr0]
            # sumcl = sumcl + (mult * (self.conc[comp][yr] - self.conc[comp][yr0])) ** 1.7
            sumcl = sumcl + np.sign(num) * np.power(mult * abs(num), 1.7)
        for comp, mult in self.pamset["br_dict"].items():
            sumbr = sumbr + mult * (conc[comp][yr] - conc[comp][yr0])
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
        yr_0 = self.years[0]
        tracer = "TROP_O3"
        # ALOG(1700.0))  !Concentration in 2010 &
        self.conc[tracer][yr] = (
            30.0
            + 6.7
            * (
                np.log(self.conc["CH4"][yr])
                - np.log(self.conc_in["CH4"][self.pamset["ref_yr"]])
            )
            + 0.17 * (self.emis["NOx"][yr] - self.emis["NOx"][self.pamset["ref_yr"]])
            + 0.0014 * (self.emis["CO"][yr] - self.emis["CO"][self.pamset["ref_yr"]])
            + 0.0042
            * (self.emis["NMVOC"][yr] - self.emis["NMVOC"][self.pamset["ref_yr"]])
        )
        # RBS101115
        # IF (yr_ix.LT.yr_2010) THEN ! Proportional to TROP_O3 build-up
        # Rewritten a bit, place to check for differences...
        value = self.conc[tracer][yr]
        value_0 = self.conc[tracer][yr_0]
        if value_0 == 30.0:
            return 0.0
        q = (value - value_0) / (30.0 - value_0) * (self.pamset["qo3"])
        return q

    def calc_aerosol_forcing(self, yr, tracer):
        """
        Calculate aerosol forcing for single aerosol

        Parameters
        ----------
        yr : int
            Year for which to calculate
        tracer : str
            Name of aerosol component to calculate for

        Returns
        -------
        float
            Forcing from aerosol in question for year in question
        """
        ref_emission_species = {
            "SO2": ["SO2", self.pamset["qdirso2"]],
            "SO4_IND": ["SO2", self.pamset["qindso2"]],
            "OC": ["OC", self.pamset["qoc"]],
            "BC": ["BC", self.pamset["qbc"]],
            "BMB_AEROS": ["BMB_AEROS_OC", self.pamset["qbmb"]],
            "NMVOC": ["NMVOC", self.pamset["qnmvoc"]],
            "NH3": ["NH3", self.pamset["qnh3"]],
            "NOx": ["NOx", self.pamset["qnox"]],
        }
        # Natural emissions
        # (after IPCC TPII on simple climate models, 1997)
        # enat = 42.0 Not used, why is this here?
        # Aerosol forcing used to be scaled to reference year
        # No just total forcing change
        # Only with emission to concentration factors differing
        # These are held in dictionary

        yr_0 = self.years[0]
        q = 0
        # Natural emissions
        # (after IPCC TPII on simple climate models, 1997)
        # enat = 42.0 Not used, why is this here?
        # Emission in reference year
        # SO2, SO4_IND, BC and OC are treated exactly the same
        # Only with emission to concentration factors differing
        # These are held in dictionary
        if self.df_gas["ALPHA"][tracer] != 0:
            q = (self.emis[tracer][yr] - self.emis[tracer][yr_0]) * self.df_gas[
                "ALPHA"
            ][tracer]
        elif tracer in ref_emission_species:
            em_change = (
                self.emis[ref_emission_species[tracer][0]][yr]
                - self.emis[ref_emission_species[tracer][0]][yr_0]
            )
            q = ref_emission_species[tracer][1] * em_change
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
          Year for which to calculate
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
        # Intialising with the combined values from CO2, N2O and CH4
        tot_forc, forc_nh, forc_sh = self.calculate_forc_three_main(yr)
        yr_0 = self.years[0]
        # Finish per tracer calculations, add per tracer to printable df and sum the total
        for tracer, forc_val_series in self.forc.items():
            if tracer in ["CO2", "N2O", "CH4"]:
                continue
            q = 0
            if tracer == "LANDUSE":
                q = rf_luc

            elif tracer in ["NMVOC", "NH3", "NOx"] or tracer[:2] in [
                "BC",
                "OC",
                "SO",
                "BM",
            ]:
                q = self.calc_aerosol_forcing(yr, tracer)
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
            forc_val_series[yr - yr_0] = q  # + FORC_PERT(yr_ix,trc_ix)
            # Calculating hemispheric forcing:
            forc_nh, forc_sh = calculate_hemispheric_forcing(
                tracer, q, forc_nh, forc_sh
            )
            tot_forc = tot_forc + q
            # print("Forcer: %s, tot_forc: %f, FN: %f, FS: %f, q: %f"%(tracer, tot_forc, forc_nh, forc_sh, q)

        # Add total forcing from precalculated:
        precalc_rf = self.precalc_dict["precalc_erf"]["TOT_vanilla"][yr]
        tot_forc, forc_nh, forc_sh = (
            tot_forc + precalc_rf,
            forc_nh + precalc_rf,
            forc_sh + precalc_rf,
        )

        # Adding forcing perturbations if they exist:
        # TODO: Deal with forcing perturbation if precalculated species is perturbed
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

    def emi2conc(self, yr, feedback_dict=None):
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
        feedback_dict : dict
            Dictionary containing feedback variables and their values
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
            self.conc["CO2"][yr] = self.carbon_cycle.co2em2conc(
                yr,
                self.emis["CO2_FF"][yr]
                + self.emis["CO2_AFOLU"][yr]
                + self.df_gas["NAT_EM"]["CO2"],
                feedback_dict=feedback_dict,
            )
            self.fill_one_row_conc(yr, avoid=["CO2"])
            return
        self.add_row_of_zeros_conc(yr)

        for tracer, value_dict in self.conc.items():
            # something something postscenario...?
            if self.df_gas["CONC_UNIT"][tracer] == "-":
                # Forcing calculated from emissions, should now only be TROP_O3
                continue
            if tracer == "CO2":
                value_dict[yr] = self.carbon_cycle.co2em2conc(
                    yr,
                    self.emis["CO2_FF"][yr]
                    + self.emis["CO2_AFOLU"][yr]
                    + self.df_gas["NAT_EM"]["CO2"],
                    feedback_dict=feedback_dict,
                )
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
            value_dict[yr] = point_conc / q + (conc_local - point_conc / q) * np.exp(-q)

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

    def write_output_to_files(self, cfg, feedback_dict_series=None, make_plot=False):
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
        feedback_dict_series : dict
            Dictionary containing feedback data for each key
        make_plot : bool
           Whether the output should be plottet or not
        """
        if "output_folder" in cfg:
            # Make os independent?
            outdir = os.path.join(os.getcwd(), cfg["output_folder"])
        else:
            outdir = os.path.join(os.getcwd(), "output")

        results_dict = self.add_results_to_dict(
            cfg, feedback_dict_series=feedback_dict_series
        )

        if "output_prefix" in cfg:
            filename_start = cfg["output_prefix"]
        else:
            filename_start = "output"

        longname_shortname_dict = {
            "emissions": "em",
            "concentrations": "conc",
            "forcing": "forc",
            "carbon cycle": "carbon",
        }
        for outtype, df in results_dict.items():

            df.to_csv(
                os.path.join(
                    outdir, f"{filename_start}_{longname_shortname_dict[outtype]}.txt"
                ),
                sep="\t",
                index=False,
                float_format="%.5e",
            )

        if make_plot:
            plot_output2("forc", results_dict["forcing"], outdir)
            plot_output2(
                "emis", results_dict["emissions"], outdir, self.df_gas["EM_UNIT"]
            )
            plot_output2(
                "conc", results_dict["concentrations"], outdir, self.df_gas["CONC_UNIT"]
            )

        if "carbon_cycle_outputs" in cfg:
            # Adding carbon cycle outputs here
            # Typically back_calculated emissions for conc_run
            # Airborne fraction
            # biosphere carbon flux
            # Ocean carbon flux
            # Yearly fluxes

            results_dict["carbon cycle"].to_csv(
                os.path.join(outdir, f"{filename_start}_carbon.txt"),
                sep="\t",
                index=False,
                float_format="%.5e",
            )

    def add_results_to_dict(self, cfg, feedback_dict_series=None):
        """
        Add results to results dictionary

        Parameters
        ----------
        cfg : dict
            Configurations to define where to put output
            files and what prefix to have for file name
            At the moment this method only needs to know
            if it's supposed to include carbon cycle outputs
        feedback_dict_series : dict
            Dictionary containing yearly feedback variables
            for each key to get feedback adjusted carbon cycle
            outputs

        Returns
        -------
        dict
            Containing run emissions, concentrations and forcings
            in the form of pandas.Dataframe with years and tracers
        """
        df_forc = pd.DataFrame(data=self.forc, index=self.years)
        df_forc["Year"] = self.years
        # Adding in precalculated values:
        df_forc = pd.concat([df_forc, self.precalc_dict["precalc_erf"]], axis=1)
        df_forc.drop(columns=["TOT_vanilla"], inplace=True)

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
        df_conc = pd.concat([df_conc, self.precalc_dict["precalc_conc"]], axis=1)
        cols = df_conc.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_conc = df_conc[cols]
        frame_order = self.df_gas.index.copy()
        frame_order = frame_order.drop(list(set(frame_order.to_list()) - set(cols)))
        frame_order = frame_order.insert(0, "Year")

        df_conc = df_conc[frame_order]

        results = {}
        results["emissions"] = df_emis
        results["concentrations"] = df_conc
        results["forcing"] = df_forc

        if "carbon_cycle_outputs" in cfg:
            results["carbon cycle"] = self.get_carbon_cycle_data(
                feedback_dict_series=feedback_dict_series
            )
        return results

    def get_carbon_cycle_data(self, feedback_dict_series=None):
        """
        Get carbon cycle data and put in dataframe for output

        Parameters
        ----------
        feedback_dict_series : dict
            Dictionary containing yearly feedback variables
            for each key to get feedback adjusted carbon cycle
            outputs

        Returns
        -------
            Pandas.DataFrame
            With carbon cycle inputs including Airborne fraction
            backcalculated emissions (in the case of concenration runs)
            Biosphere carbon flux and ocean carbon flux
        """
        conc_series = np.array([v for k, v in self.conc["CO2"].items()])
        if self.pamset["conc_run"]:
            df_carbon = self.carbon_cycle.get_carbon_cycle_output(
                self.years,
                conc_series=conc_series,
                conc_run=self.pamset["conc_run"],
                feedback_dict_series=feedback_dict_series,
            )
        else:
            em_series = (
                self.emis["CO2_FF"][self.years].values
                + self.emis["CO2_AFOLU"][self.years].values
            )
            df_carbon = self.carbon_cycle.get_carbon_cycle_output(
                self.years,
                feedback_dict_series=feedback_dict_series,
                conc_series=conc_series,
            )
            if df_carbon is None:
                df_carbon = pd.DataFrame(
                    data={
                        "Airborne fraction CO2": calculate_airborne_fraction(
                            em_series, conc_series
                        )
                    },
                    index=self.years,
                )
            else:
                df_carbon["Airborne fraction CO2"] = calculate_airborne_fraction(
                    em_series, conc_series
                )
        return df_carbon

    def get_feedback_list(self):
        """
        Retrieve the list of feedback variable names from the carbon cycle.

        Returns
        -------
        list
            A list containing the feedback variable names that the carbon cycle model uses.
        """
        return self.carbon_cycle.get_feedback_list()
