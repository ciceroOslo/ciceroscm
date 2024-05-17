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

# Conversion factor ppm CO2 -> kg
PPM_CO2_TO_PG_C = 2.123


def xco2_poly_to_solve(z_co2, constant=0, mixed_carbon=0):
    """
    Define function to invert to get zco2 from xco2y assuming 0th step
    """
    return (
        (1.3021 + 2 * mixed_carbon / GE_COEFF / PPMKG_TO_UMOL_PER_VOL) * z_co2
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

    Parameters
    ----------
    pamset : dict
        parameter set dictionary which might contain the value "missing"
        for the rs_function and rb_function because of the way we deal
        with expected values in the parameter set. In this case
        we need to take them out of the parameterset so they can be run
        with defaults

    Returns
    -------
    dict
        Updated version of the parameterset with "missing" values deleted
    """
    for key, value in pamset.items():
        if value == "missing":
            del pamset[key]
    return pamset


def calculate_airborne_fraction(em_timeseries, conc_timeseries):
    """
    Calculate Airborne Fraction of CO2 from emissions timeseries

    Parameters
    ----------
    em_timeseries: np.ndarray
        Emissions timeseries, either inputs or backcalculated for concentration
        run
    conc_timeseries : np.ndarray
        Concentrations timeseries. Should be the same length as the emissions
        timeseries

    Returns
    -------
    np.ndarray
        Airborne fraction calculated from the em_timeseries and conc_timeseries
    """
    airborne_fraction = (
        (conc_timeseries - 278.0) / np.cumsum(em_timeseries) * PPM_CO2_TO_PG_C
    )
    return airborne_fraction


class CarbonCycleModel:
    """
    Class to handle carbon cycle calculations
    """

    def __init__(self, pamset):
        """
        Initialise Carbon cycle model

        Parameters
        ----------
            pamset : dict
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

    def reset_co2_hold(self, beta_f=0.287, mixed_carbon=75.0):
        """
        Reset values of CO2_hold for new run

        This method is mainly called to do a new run with the same cscm instance,
        in which case you need to reset hold values, and be able to update
        parameter values for the carbon cycle free parameter beta_f
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
        self.pamset["mixed_carbon"] = mixed_carbon

    def _set_co2_hold(
        self, xco2=278.0, yco2=0.0, emco2_prev=0.0, ss1=0.0, sums=0
    ):  # pylint: disable=too-many-arguments
        """
        Reset the CO2 hold scalar values,

        Use this to rerun from same state as before
        in year not zero. Should only be used for back-calculations

        Parameters
        ----------
        xco2 : float
            CO2 concentration to set, default is 278.0 which is the start value
        yco2 : float
            yco2 value, default is 0.0 which is the start value
        emco2_prev : float
            emissions in previous timestep, default is 0.0
        ss1 : float
            ss1 value, default is 0.0
        sums : float
            sums of ocean uptake inorganic carbon, default is 0.0
        """
        self.co2_hold["yCO2"] = yco2
        self.co2_hold["xCO2"] = xco2
        self.co2_hold["emCO2_prev"] = emco2_prev
        self.co2_hold["ss1"] = ss1
        self.co2_hold["sums"] = sums

    def _get_co2_hold_values(self):
        """
        Get co2_hold scalar values as dictionary

        These can be used to reset the values after, this is useful for
        the back calculation which needs to run the same timestep over
        and over. This is why the method is private as it is just meant
        to be called from the back-calculations to get values that can then
        be sent to the _set_co2 hold function
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
            em_co2 = (em_co2_common - ffer) / PPM_CO2_TO_PG_C

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
                / self.pamset["mixed_carbon"]
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
        Get the biospheric fertilisation (ffer) time series

        For emissions runs this yields appropriate amounts only for years for which the
        model has been run and zeros otherwise. If a concentrations series is used
        this yields a back calculated estimate of the carbon pool the model
        estimates given these atmospheric concentrations.

        Parameters
        ----------
        conc_run : bool
            Whether calculation is being done for a conc_run
            (i.e. a run that has CO2_hold values precalculated, or whether back
            calculations are needed or co2_conc_series need to be used for calculations)
        co2_conc_series : np.ndarray
            Time series of concentrations to calculate ffer from

        Returns
        -------
        np.ndarray
            Biospheric fertilisation factor timeseries
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

        Parameters
        ----------
        conc_run : bool
            Whether calculation is being done for a conc_run
            (i.e. a run that has CO2_hold values precalculated, or whether back
            calculations are needed or co2_conc_series need to be used for calculations)
        co2_conc_series : np.ndarray
            Time series of concentrations to calculate ffer from

        Returns
        -------
        np.ndarray
            Timeseries of the added carbon content to the biosphere carbon pool

        # TODO: Make this be flux rather than cumulative
        """
        ffer = self._get_ffer_timeseries(conc_run, co2_conc_series)
        biosphere_carbon_pool = (
            np.array(
                [
                    np.sum(ffer[: self.pamset["idtm"] * (yrix + 1)])
                    for yrix in range(self.pamset["years_tot"])
                ]
            )
            / self.pamset["idtm"]
        )
        return biosphere_carbon_pool

    def get_ocean_carbon_pool_content(self, conc_run=False, co2_conc_series=None):
        """
        Get ocean carbon pool content

        If this is called from a concentrations either it is
        assumed that emissions are already back calculated
        Otherwhise you can send a concentration series to do
        the back calculations with to set all the values in
        the sCO2 array

        Returns
        -------
        np.ndarray
            Timeseries of the added carbon content to the ocean carbon pool

        TODO: Have this be flux rather than cumulative
        """
        # ocean_carbon_pool = np.array([
        #    np.sum(self.co2_hold["sCO2"][: self.pamset["idtm"] * (yrix+1)])
        #    for yrix in range(self.pamset["years_tot"])
        # ])
        if conc_run and co2_conc_series is not None:
            self.back_calculate_emissions(co2_conc_series)
        ocean_carbon_pool = (
            np.array(
                [
                    np.sum(self.co2_hold["sCO2"][: self.pamset["idtm"] * (yrix + 1)])
                    for yrix in range(self.pamset["years_tot"])
                ]
            )
            * PPM_CO2_TO_PG_C
            / self.pamset["idtm"]
            * GE_COEFF
            * OCEAN_AREA
        )
        return ocean_carbon_pool

    def back_calculate_emissions(self, co2_conc_series):
        """
        Back calculate emissions from conc_run

        co2_conc_series is assumed to be the series of concentrations
        from the year 0

        Parameters
        ----------
        co2_conc_series : np.ndarray
            Timeseries of co2 concentrations for which to back
            calculate emissions
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
                self._guess_emissions_iteration(
                    co2_conc, prev_co2_conc, yrix=i, ffer=ffer_here
                )
                * PPM_CO2_TO_PG_C
                + ffer_here
            )
            prev_co2_conc = co2_conc
        return em_series

    def simplified_em_backward(self, co2_conc_now, co2_conc_zero):
        """
        Simplified _guess solution to find emissions

        Based on algebraic solution for single timestep which is only
        valid for the very first timestep

        Parameters
        ----------
        co2_conc_now : float
            Value of CO2 concentration in timestep resulting after
            the emissions you would like to find are applied
        co2_conc_zero : float
            Value of CO2 concentration in timestep before the step
            for which you want to find the concentrations

        Returns
        -------
        float
            Emissions from the simplified algebraic approach
        # TODO : Should this be private?
        """
        co2_diff = co2_conc_zero - co2_conc_now
        poly_of_z = partial(xco2_poly_to_solve, constant=co2_diff)
        z_solve = optimize.fsolve(poly_of_z, 0)
        cc1 = OCEAN_AREA * GE_COEFF / (1 + OCEAN_AREA * GE_COEFF / 2.0)
        return (
            z_solve
            * 2
            * self.pamset["mixed_carbon"]
            / PPMKG_TO_UMOL_PER_VOL
            * OCEAN_AREA
            / cc1
        )

    def _guess_emissions_iteration(
        self, co2_conc_now, co2_conc_zero, yrix=0, rtol=1e-7, maxit=100, ffer=None
    ):  # pylint: disable=too-many-arguments
        """
        Iterate to get right emissions for a single year

        Make sure the yrix is where you are at

        Parameters
        ----------
        co2_conc_now : float
            Value of CO2 concentration in timestep resulting after
            the emissions you would like to find are applied
        co2_conc_zero : float
            Value of CO2 concentration in timestep before the step
            for which you want to find the concentrations
        yrix : int
            Index of year to do back calculation for, should be the
            year number starting with 0 of the total of years for which
            the back calculated emissions in total for which this calculation
            should be made
        rtol : float
            Relative tolerance in accuracy between concentrations change from
            back calculated emissions set and input concentrations
        maxit : int
            Maximum number of iterations todo before cutting. This is a
            safety switch to make sure we don't do infinite looping if
            the solution doesn't converge
        ffer : np.ndarray
            biospheric fertilisation for the given concentrations change

        Returns
        -------
        float
            Back calculated emissions that yields a concentrations change
            from the conc_co2_zero to conc_co2_now which

        """
        if ffer is None:
            ffer = self._get_ffer_timeseries([co2_conc_zero, co2_conc_now])[
                yrix * self.pamset["idtm"]
            ]
        min_guess = self.simplified_em_backward(co2_conc_now / 2, co2_conc_zero)
        max_guess = self.simplified_em_backward(co2_conc_now * 2, co2_conc_zero)
        guess = self.simplified_em_backward(co2_conc_now, co2_conc_zero)
        hold_dict = self._get_co2_hold_values()
        estimated_conc = self.co2em2conc(
            self.pamset["nystart"] + yrix, PPM_CO2_TO_PG_C * (guess) + ffer
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
                self.pamset["nystart"] + yrix, PPM_CO2_TO_PG_C * (guess) + ffer
            )
            iteration = iteration + 1
        return guess
