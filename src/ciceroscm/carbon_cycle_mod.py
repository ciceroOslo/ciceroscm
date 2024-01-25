"""
Module to handle carbon cycle from CO2 emissions to concentrations
"""

import numpy as np

from ._utils import cut_and_check_pamset


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

    def co2em2conc(self, yr, em_co2_common):  # pylint: disable=too-many-locals
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
        # Area of the ocean (m^2)
        ocean_area = 3.62e14

        # Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
        coeff = 1.0 / (ocean_area * 9.06)

        # TIMESTEP (YR)
        dt = 1.0 / self.pamset["idtm"]

        # Conversion factor ppm/kg --> umol*/m3
        conv_factor = 1.722e17

        # USING MIXED LAYER DEPTH = 75 metres # TODO: Should this be synced to what's in the UDM?
        mixed_layer_depth = 75.0

        cc1 = dt * ocean_area * coeff / (1 + dt * ocean_area * coeff / 2.0)
        yr_ix = yr - self.pamset["nystart"]
        # Monthloop:
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
        return self.co2_hold["xCO2"]
