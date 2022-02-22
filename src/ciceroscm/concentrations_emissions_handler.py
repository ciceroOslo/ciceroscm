"""
Module to handle emission, concentrations and converstion to forcing
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
    """
    df_gas = pd.read_csv(filename, delim_whitespace=True, index_col=0)
    df_gas.rename(
        columns={"TAU1(YEARS)": "TAU1", "NATURAL_EMISSIONS": "NAT_EM"}, inplace=True
    )
    return df_gas


def _rs_function(it, idtm=24):
    """
    Calculate pulse response function for mixed layer
    time is the year index*idtm + i, i.e. the month number
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
    time is the year index*idtm + i, i.e. the month number
    """
    time = it / idtm
    biotic_decay = (
        0.70211 * np.exp(-0.35 * time)
        + 13.4141e-3 * np.exp(-time / 20.0)
        - 0.71846 * np.exp(-55 * time / 120.0)
        + 2.9323e-3 * np.exp(-time / 100.0)
    )
    return biotic_decay


def read_natural_emissions(filename, component, startyear=1750, endyear=2500):
    """
    Read in single column natural emissions file
    """
    df_natemis = pd.read_csv(filename, header=None, names=[component], index_col=False)
    df_natemis["year"] = np.arange(startyear, endyear + 1)
    df_natemis = df_natemis.set_index("year")
    return df_natemis


def check_pamset(pamset):
    """
    Check that parameterset has necessary values for run
    Otherwise set to default values
    """
    required = {
        "lamb": 0.8,
        "qbmb": 0.03,
        "qo3": 0.4,
        "qdirso2": -0.457,
        "qindso2": -0.514,
        "qbc": 0.200,
        "qoc": -0.103,
        "ref_yr": 2010,
        "idtm": 24,
    }

    pamset = check_numeric_pamset(required, pamset)
    if "lifetime_mode" not in pamset:
        pamset["lifetime_mode"] = "TAR"

    return pamset


def read_inputfile(input_file):
    """
    Read input from emissions or concentrations file
    """
    df_input = pd.read_csv(
        input_file, delim_whitespace=True, index_col=0, skiprows=[1, 2, 3]
    )
    return df_input


def calculate_hemispheric_forcing(tracer, q, forc_nh, forc_sh):
    """
    Calculate hemispheric forcing per tracer
    """
    if tracer in ("SO2", "SO4_IND"):
        forc_nh = forc_nh + q * 1.6
        forc_sh = forc_sh + q * 0.4
    elif tracer == "TROP_O3":
        # 1.29+0.74 != 2, update/check
        forc_nh = forc_nh + q * 1.29  # 1.3
        forc_sh = forc_sh + q * 0.74  # 0.7
    else:
        forc_nh = forc_nh + q
        forc_sh = forc_sh + q
    return forc_nh, forc_sh


class ConcentrationsEmissionsHandler:
    """
    Class to handle concentrations
    and emissions input for ciceroscm
    """

    # pylint: disable=too-many-instance-attributes
    # Need a few more here or else consider breaking
    # up emissions in it's own class or
    # CO2-handling in its own class

    def __init__(
        self, cfg, pamset,
    ):
        self.df_gas = read_components(cfg["gaspamfile"])
        self.conc = {}
        self.forc = {}
        if "concentrations_file" in cfg:
            self.conc_in = read_inputfile(cfg["concentrations_file"])
        if "emissions_file" in cfg:
            self.emis = read_inputfile(cfg["emissions_file"])
            for tracer in self.df_gas.index:
                if tracer != "CO2.1":
                    self.conc[tracer] = {}
                    self.forc[tracer] = []
            self.forc["Total_forcing"] = []
            self.emis.rename(
                columns={"CO2": "CO2_FF", "CO2.1": "CO2_AFOLU"}, inplace=True
            )
        self.nat_emis_ch4 = read_natural_emissions(cfg["nat_ch4_file"], "CH4")
        self.nat_emis_n2o = read_natural_emissions(cfg["nat_n2o_file"], "N2O")
        self.pamset = check_pamset(pamset)
        if "conc_run" in cfg:
            self.pamset["conc_run"] = cfg["conc_run"]
        else:
            self.pamset["conc_run"] = False
        self.years = np.arange(self.pamset["nystart"], self.pamset["nyend"] + 1)

        self.co2_hold = {
            "yCO2": 0.0,
            "xCO2": 278.0,
            "sCO2": np.zeros(self.pamset["idtm"] * len(self.years)),
            "emCO2_prev": 0.0,
            "dfnpp": np.zeros(self.pamset["idtm"] * len(self.years)),
            "ss1": 0.0,
            "sums": 0.0,
        }
        self.co2_hold["dfnpp"]
        self.r_functions = np.empty(
            (2, self.pamset["idtm"] * 351)
        )  # if speedup, get this to reflect number of years
        self.r_functions[0, :] = [
            _rs_function(it, self.pamset["idtm"])
            for it in range(self.pamset["idtm"] * len(self.years))
        ]
        self.r_functions[1, :] = [
            _rb_function(it, self.pamset["idtm"])
            for it in range(self.pamset["idtm"] * len(self.years))
        ]

    def calculate_strat_quantities(self, yr):
        """
        Calculate sumcl and sumbr in stratosphere
        sums of effects on stratospheric ozone from
        chlorinated gases and halons
        """
        if yr < self.years[0]:
            return 0.0, 0.0
        chlor_dict = {
            "CFC-11": 3,
            "CFC-12": 2,
            "CFC-113": 3,
            "CFC-114": 2,
            "CFC-115": 1,
            "CCl4": 4,
            "CH3CCl3": 3,
            "HCFC-22": 1,
            "HCFC-141b": 2,
            "HCFC-123": 2,
            "HCFC-142b": 1,
        }
        # More Halons?
        brom_dict = {"H-1211": 1, "H-1301": 1}
        sumcl = 0
        sumbr = 0
        for comp, mult in chlor_dict.items():
            sumcl = sumcl + (mult * self.conc[comp][yr]) ** 1.7
        for comp, mult in brom_dict.items():
            sumbr = sumbr + (mult * self.conc[comp][yr])
        return sumcl, sumbr

    def calculate_forc_three_main(self, yr):
        """
        Calculate forcings for components CO2, N2O and CH4
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
        c0_n2o = self.conc["N2O"][self.years[0]]
        c_n2o = self.conc["N2O"][yr]
        c0_co2 = self.conc["CO2"][self.years[0]]
        c_co2 = self.conc["CO2"][yr]
        c0_ch4 = self.conc["CH4"][self.years[0]]
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
        q_ch4 = 1.0 / 1.14 * q_ch4  # + FORC_PERT(yr_ix,trc_ix))
        self.forc["CO2"].append(q_co2)
        self.forc["CH4"].append(q_ch4)
        self.forc["N2O"].append(q_n2o)
        tot_forc = q_co2 + q_n2o + q_ch4
        forc_nh, forc_sh = calculate_hemispheric_forcing(
            "three_main", tot_forc, 0.0, 0.0
        )

        return tot_forc, forc_nh, forc_sh

    def tropospheric_ozone_forcing(self, yr):
        """
        Calculate tropospheric ozone forcing
        """
        emstart = self.pamset["emstart"]
        yr_emstart = emstart - self.years[0]
        yr_ix = yr - self.years[0]
        tracer = "TROP_O3"
        if self.pamset["conc_run"] or yr_ix < yr_emstart:
            # Uses change in CO2_FF emissions
            q = (
                (self.emis["CO2_FF"][yr] - self.emis["CO2_FF"][self.years[0]])
                / (
                    self.emis["CO2_FF"][self.pamset["ref_yr"]]
                    - self.emis["CO2_FF"][self.years[0]]
                )
                * self.pamset["qo3"]
            )

        else:
            # ALOG(1700.0))  !Concentration in 2010 &
            self.conc[tracer][yr] = (
                30.0
                + 6.7 * (np.log(self.conc["CH4"][yr]) - np.log(1832.0))
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

    def conc2forc(self, yr, rf_luc, rf_sun):
        """
        Calculate forcing from concentrations
        """
        ref_emission_species = {
            "SO2": ["SO2", self.pamset["qdirso2"]],
            "SO4_IND": ["SO2", self.pamset["qindso2"]],
            "OC": ["OC", self.pamset["qoc"]],
            "BC": ["BC", self.pamset["qbc"]],
        }
        # Intialising with the combined values from CO2, N2O and CH4
        tot_forc, forc_nh, forc_sh = self.calculate_forc_three_main(yr)
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
                erefyr = self.emis[ref_emission_species[tracer][0]][
                    self.pamset["ref_yr"]
                ]
                if erefyr != 0.0:  # pylint: disable=compare-to-zero
                    frac_em = self.emis[ref_emission_species[tracer][0]][yr] / erefyr
                    q = ref_emission_species[tracer][1] * frac_em

            elif (
                tracer in self.df_gas.index
                and self.df_gas["ALPHA"][tracer] != 0  # pylint: disable=compare-to-zero
            ):
                q = (
                    (self.conc[tracer][yr] - self.conc[tracer][self.years[0]])
                    * self.df_gas["ALPHA"][tracer]
                )  # +forc_pert
            elif tracer == "TROP_O3":
                q = self.tropospheric_ozone_forcing(yr)
            elif tracer == "STRAT_O3":
                sumcl, sumbr = self.calculate_strat_quantities(yr - 3)
                q = -(0.287737 * (0.000552 * (sumcl) + 3.048 * sumbr)) / 1000.0
            elif tracer == "STRAT_H2O":
                q = (
                    0.15 * 1.14 * self.forc["CH4"][yr - self.years[0]]
                )  # + FORC_PERT(yr_ix,trc_ix)
            elif tracer == "OTHER":
                # Possible with forcing perturbations for other
                # components such as contrails, cirrus etc...
                pass
            self.forc[tracer].append(q)  # + FORC_PERT(yr_ix,trc_ix)
            # Calculating hemispheric forcing:
            forc_nh, forc_sh = calculate_hemispheric_forcing(
                tracer, q, forc_nh, forc_sh
            )
            tot_forc = tot_forc + q
            # print("Forcer: %s, tot_forc: %f, FN: %f, FS: %f, q: %f"%(tracer, tot_forc, forc_nh, forc_sh, q))
        # Adding solar forcing
        # tot_forc = tot_forc + rf_sun
        self.forc["Total_forcing"].append(tot_forc)
        forc_nh = forc_nh + rf_sun
        forc_sh = forc_sh + rf_sun
        # print("yr: %d, tot_forc: %f, FN: %f, FS: %f "%(yr, tot_forc, forc_nh, forc_sh))

        return tot_forc, forc_nh, forc_sh

    def emi2conc(self, yr):
        """
        Calculate concentrations from emissions
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
                self.df_gas["NAT_EM"][tracer] = self.nat_emis_ch4["CH4"][yr]
                q = self.methane_lifetime(q, conc_local, yr)
            if tracer == "N2O":
                self.df_gas["NAT_EM"][tracer] = self.nat_emis_n2o["N2O"][yr]

            emis = self.emis[tracer][yr]
            emis = (
                emis + self.df_gas["NAT_EM"][tracer]
            )  # natural emissions, from gasspamfile
            point_conc = emis / self.df_gas["BETA"][tracer]
            # Rewrote this quite a bit from an original loop,
            # but I think it is mathematically equivalent
            self.conc[tracer][yr] = point_conc / q + (
                conc_local - point_conc / q
            ) * np.exp(-q)

    def methane_lifetime(self, q, conc_local, yr):
        """
        Calculate methane concentrations from emissions
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
        elif self.pamset["lifetime_mode"] == "CONSTANT":
            q = 1.0 / 12.0
        else:
            q = q * (((conc_local / 1700.0)) ** (ch4_wigley_exp))

        q = q + 1.0 / self.df_gas["TAU2"]["CH4"] + 1.0 / self.df_gas["TAU3"]["CH4"]

        return q

    def co2em2conc(self, yr):  # pylint: disable=too-many-locals
        """
        Calculate co2 concentrations from emissions
        """
        # Fertilisation factor
        beta_f = 0.287

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
                    60 * beta_f * np.log(self.co2_hold["xCO2"] / 278.0)
                )
            if it > 0:
                for j in range(1, it):
                    sumf = (
                        sumf
                        + self.co2_hold["dfnpp"][j] * self.r_functions[1, it - 1 - j]
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
            sumz = 0.0
            for j in range(it - 1):
                sumz = sumz + self.co2_hold["sCO2"][j] * self.r_functions[0, it - 1 - j]
            z_co2 = (
                conv_factor
                * coeff
                * dt
                / mixed_layer_depth
                * (sumz + 0.5 * self.co2_hold["sCO2"][it])
            )
            self.co2_hold["yCO2"] = (
                1.3021 * z_co2
                + 3.7929e-3 * (z_co2 ** 2)
                + 9.1193e-6 * (z_co2 ** 3)
                + 1.488e-8 * (z_co2 ** 4)
                + 1.2425e-10 * (z_co2 ** 5)
            )
            self.co2_hold["xCO2"] = (
                self.co2_hold["sCO2"][it] + self.co2_hold["yCO2"] + 278.0
            )
            # print("it: %d, emCO2: %e, sCO2: %e, zCO2: %e, yCO2: %e, xCO2: %e, ss1: %e, ss2: %e, dnfpp:%e"%(it, em_co2, self.co2_hold["sCO2"][it], z_co2, self.co2_hold["yCO2"], self.co2_hold["xCO2"], self.co2_hold["ss1"], ss2, self.co2_hold["dfnpp"][it]))
        self.conc["CO2"][yr] = self.co2_hold["xCO2"]

    def fill_one_row_conc(
        self, yr, avoid=[]
    ):  # pylint: disable=dangerous-default-value
        """
        Fill in one row of concentrations in conc_dict
        """
        for tracer, value_dict in self.conc.items():
            if tracer in avoid:
                continue
            if tracer in self.conc_in:
                value_dict[yr] = self.conc_in[tracer][yr]
            else:
                value_dict[yr] = 0

    def add_row_of_zeros_conc(self, yr):
        """
        Fill in one row of concentrations in conc_dict
        """
        for value_dict in self.conc.values():
            value_dict[yr] = 0

    def write_output_to_files(self, cfg):
        """
        Write results to files after run
        """
        if "output_prefix" in cfg:
            # Make os independent?
            outdir = os.path.join(os.getcwd(), cfg["output_prefix"])
        else:
            outdir = os.path.join(os.getcwd(), "output")

        df_forc = pd.DataFrame(data=self.forc, index=self.years)
        df_forc["Year"] = self.years

        cols = df_forc.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df_forc = df_forc[cols]
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
        df_forc.to_csv(
            os.path.join(outdir, "output_forc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
        df_conc.to_csv(
            os.path.join(outdir, "output_conc.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )

        df_emis.to_csv(
            os.path.join(outdir, "output_em.txt"),
            sep="\t",
            index=False,
            float_format="%.5e",
        )
