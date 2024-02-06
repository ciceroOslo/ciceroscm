"""
Module to handle carbon cycle from CO2 emissions to concentrations
"""

from functools import partial

import numpy as np
from scipy import optimize

from ._utils import cut_and_check_pamset

# Area of the ocean (m^2)
OCEAN_AREA = 3.62e14

# Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
GE_COEFF = 1.0 / (OCEAN_AREA * 9.06)


# Conversion factor ppm/kg --> umol*/m3
PPMKG_TO_UMOL_PER_VOL = 1.722e17

# USING MIXED LAYER DEPTH = 75 metres # TODO: Should this be synced to what's in the UDM?
MIXED_LAYER_DEPTH = 75.0


def xco2_poly_to_solve(z_co2, constant=0):
    """
    Define function to invert to get zco2 from xco2y assuming 0th step
    """
    return (
        (1.3021 + 2 * MIXED_LAYER_DEPTH / GE_COEFF / PPMKG_TO_UMOL_PER_VOL) * z_co2
        + 3.7929e-3 * (z_co2**2)
        + 9.1193e-6 * (z_co2**3)
        + 1.488e-8 * (z_co2**4)
        + 1.2425e-10 * (z_co2**5)
        + constant
    )


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


def take_out_missing(pamset):
    """
    Take out values that are missing from pamset

    Needed to take care of rs_function and rb_function
    """
    for key, value in pamset.items():
        if value == "missing":
            del pamset[key]
    return pamset


class CarbonCycleModel:
    """
    Class to handle carbon cycle calculations
    """

    def __init__(self, pamset):
        """
        Initialise Carbon cycle model
        """
        pamset = take_out_missing(pamset.copy())
        self.pamset = cut_and_check_pamset(
            {"idtm": 24, "nystart": 1750},
            pamset,
            used={"rs_function": _rs_function, "rb_function": _rb_function},
        )
        self.pamset["years_tot"] = pamset["nyend"] - pamset["nystart"] + 1
        self.reset_co2_hold()
        self.precalc_r_functions()

    def reset_co2_hold(self, beta_f=0.287):
        """
        Reset values of CO2_hold for new run
        """
        self.co2_hold = {
            "yCO2": 0.0,
            "xCO2": 278.0,
            "sCO2": np.zeros(self.pamset["idtm"] * self.pamset["years_tot"]),
            "emCO2_prev": 0.0,
            "dfnpp": np.zeros(self.pamset["idtm"] * self.pamset["years_tot"]),
            "ss1": 0.0,
            "sums": 0.0,
        }
        self.pamset["beta_f"] = beta_f

    def _set_co2_hold(
        self, xco2=278.0, yco2=0.0, emco2_prev=0.0, ss1=0.0, sums=0
    ):  # pylint: disable=too-many-arguments
        """
        Reset the CO2 hold scalar values,

        Use this to rerun from same state as before
        in year not zero. Should only be used for back-calculations
        """
        self.co2_hold["yCO2"] = yco2
        self.co2_hold["xCO2"] = xco2
        self.co2_hold["emCO2_prev"] = emco2_prev
        self.co2_hold["ss1"] = ss1
        self.co2_hold["sums"] = sums

    def _get_co2_hold_values(self):
        """
        Get co2_hold scalar values as dictionary
        """
        scalar_dict = {
            "yco2": self.co2_hold["yCO2"],
            "xco2": self.co2_hold["xCO2"],
            "emco2_prev": self.co2_hold["emCO2_prev"],
            "ss1": self.co2_hold["ss1"],
            "sums": self.co2_hold["sums"],
        }
        return scalar_dict

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
        self.r_functions = np.empty(
            (2, self.pamset["idtm"] * self.pamset["years_tot"])
        )  # if speedup, get this to reflect number of years
        if "rs_function" not in self.pamset:
            self.pamset["rs_function"] = _rs_function
        if "rb_function" not in self.pamset:
            self.pamset["rb_function"] = _rb_function
        self.r_functions[0, :] = [
            self.pamset["rs_function"](it, self.pamset["idtm"])
            for it in range(self.pamset["idtm"] * self.pamset["years_tot"])
        ]
        self.r_functions[1, :] = [
            self.pamset["rb_function"](it, self.pamset["idtm"])
            for it in range(self.pamset["idtm"] * self.pamset["years_tot"])
        ]

    def co2em2conc(self, yr, em_co2_common):
        """
        Calculate co2 concentrations from emissions

        Method to calculate co2 concentrations from emissions
        Implementing a rudimentary carbon cycle which loops over
        idtm (usually 24) timesteps a year

        Parameters
        ----------
        yr : int
          Year for which to calculate
        em_co2_common : float
             Sum of CO2 emissions from fossil fuels, land use change and natural emissions
             for the year in question

        Returns
        -------
        float
             CO2 concetrations for year in question
        """
        # TIMESTEP (YR)
        dt = 1.0 / self.pamset["idtm"]

        cc1 = dt * OCEAN_AREA * GE_COEFF / (1 + dt * OCEAN_AREA * GE_COEFF / 2.0)
        yr_ix = yr - self.pamset["nystart"]
        # Monthloop:
        for i in range(self.pamset["idtm"]):
            it = yr_ix * self.pamset["idtm"] + i
            sumf = 0.0

            # Net emissions, including biogenic fertilization effects
            if it > 0:
                # Net primary production in timestep
                self.co2_hold["dfnpp"][it] = (
                    60 * self.pamset["beta_f"] * np.log(self.co2_hold["xCO2"] / 278.0)
                )
                # Decay from previous primary production
                sumf = float(
                    np.dot(
                        self.co2_hold["dfnpp"][1:it],
                        np.flip(self.r_functions[1, : it - 1]),
                    )
                )
            # Total biospheric sink:
            ffer = self.co2_hold["dfnpp"][it] - dt * sumf

            # Sum anthropogenic and biospheric and PPMKG_TO_UMOL_PER_VOLert gC/yr --> ppm/yr
            em_co2 = (em_co2_common - ffer) / 2.123

            if it == 0:  # pylint: disable=compare-to-zero
                self.co2_hold["ss1"] = 0.5 * em_co2 / (OCEAN_AREA * GE_COEFF)
                ss2 = self.co2_hold["ss1"]
                self.co2_hold["sums"] = 0.0
            else:
                ss2 = 0.5 * em_co2 / (OCEAN_AREA * GE_COEFF) - self.co2_hold["yCO2"] / (
                    dt * OCEAN_AREA * GE_COEFF
                )
                self.co2_hold["sums"] = (
                    self.co2_hold["sums"]
                    + self.co2_hold["emCO2_prev"] / (OCEAN_AREA * GE_COEFF)
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
                PPMKG_TO_UMOL_PER_VOL
                * GE_COEFF
                * dt
                / MIXED_LAYER_DEPTH
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
        return self.co2_hold["xCO2"]

    def _get_ffer_timeseries(self, conc_run=False, co2_conc_series=None):
        """
        Get the ffer time series

        For emissions runs this yields appropriate amounts only for years for which the
        model has been run and zeros otherwise. If a concentrations series is used
        this yields a back calculated estimate of the carbon pool the model
        estimates given these atmospheric concentrations.
        """
        dt = 1.0 / self.pamset["idtm"]

        if conc_run and co2_conc_series is not None:
            timesteps = len(co2_conc_series) * self.pamset["idtm"]
            dfnpp = np.repeat(
                [
                    60 * self.pamset["beta_f"] * np.log(co2_conc / 278.0)
                    for co2_conc in co2_conc_series
                ],
                self.pamset["idtm"],
            )
        else:
            timesteps = self.pamset["idtm"] * self.pamset["years_tot"]
            dfnpp = self.co2_hold["dfnpp"]
        sumf = np.zeros(timesteps)
        print(dfnpp)
        print(len(dfnpp))
        print(timesteps)
        print(len(self.r_functions[1, :]))
        sumf[1:] = [
            float(
                np.dot(
                    dfnpp[1:it],
                    np.flip(self.r_functions[1, : it - 1]),
                )
            )
            for it in range(1, timesteps)
        ]
        ffer = dfnpp - dt * sumf
        return ffer

    def get_biosphere_carbon_pool_content(self, conc_run=False, co2_conc_series=None):
        """
        Get the carbon amount contained in the biosphere as timeseries over years

        For emissions runs this yields appropriate amounts only for years for which the
        model has been run and zeros otherwise. If a concentrations series is used
        this yields a back calculated estimate of the carbon pool the model
        estimates given these atmospheric concentrations.
        """
        ffer = self._get_ffer_timeseries(conc_run, co2_conc_series)
        biosphere_carbon_pool = [
            np.sum(ffer[: self.pamset["idtm"] * yrix])
            for yrix in range(self.pamset["years_tot"])
        ]
        return biosphere_carbon_pool

    def back_calculate_emissions(self, co2_conc_series):
        """
        Back calculate emissions from conc_run

        co2_conc_series is assumed to be the series of concentrations
        from the year 0
        """
        # Calculating fertilisation factor for all the time steps:
        # ffer = _get_ffer_timeseries(conc_run, co2_conc_series, conc_run)
        # TODO implement
        prev_co2_conc = 278.0
        em_series = np.zeros(len(co2_conc_series))
        ffer = self._get_ffer_timeseries(conc_run=True, co2_conc_series=co2_conc_series)
        for i, co2_conc in enumerate(co2_conc_series):
            ffer_here = ffer[i * self.pamset["idtm"]]
            em_series[i] = (
                self.guess_emissions_iteration(
                    co2_conc, prev_co2_conc, yrix=i, ffer=ffer_here
                )
                * 2.123
                + ffer_here
            )
            prev_co2_conc = co2_conc
        return em_series

    def simplified_em_backward(self, co2_conc_now, co2_conc_zero):
        """
        Simplified guess solution to find emissions
        """
        co2_diff = co2_conc_zero - co2_conc_now
        poly_of_z = partial(xco2_poly_to_solve, constant=co2_diff)
        z_solve = optimize.fsolve(poly_of_z, 0)
        cc1 = OCEAN_AREA * GE_COEFF / (1 + OCEAN_AREA * GE_COEFF / 2.0)
        return (
            z_solve * 2 * MIXED_LAYER_DEPTH / PPMKG_TO_UMOL_PER_VOL * OCEAN_AREA / cc1
        )

    def guess_emissions_iteration(
        self, co2_conc_now, co2_conc_zero, yrix=0, rtol=1e-7, maxit=100, ffer=None
    ):  # pylint: disable=too-many-arguments
        """
        Iterate to get right emissions

        Make sure the yrix is where you are at
        """
        if ffer is None:
            ffer = self._get_ffer_timeseries([co2_conc_zero, co2_conc_now])[
                yrix * self.pamset["idtm"]
            ]
        min_guess = self.simplified_em_backward(co2_conc_now / 2, co2_conc_zero)
        max_guess = self.simplified_em_backward(co2_conc_now * 2, co2_conc_zero)
        guess = self.simplified_em_backward(co2_conc_now, co2_conc_zero)
        hold_dict = self._get_co2_hold_values()
        # print(yrix)
        estimated_conc = self.co2em2conc(
            self.pamset["nystart"] + yrix, 2.123 * (guess) + ffer
        )
        iteration = 0
        while (
            iteration < maxit
            and np.abs(co2_conc_now - estimated_conc) / co2_conc_now > rtol
        ):
            if estimated_conc > co2_conc_now:
                max_guess = guess
                guess = (guess + min_guess) / 2

            else:
                min_guess = guess
                guess = (guess + max_guess) / 2
            self._set_co2_hold(**hold_dict)
            estimated_conc = self.co2em2conc(
                self.pamset["nystart"] + yrix, 2.123 * (guess) + ffer
            )
            iteration = iteration + 1
        # print(yrix)
        # print(guess)
        return guess
