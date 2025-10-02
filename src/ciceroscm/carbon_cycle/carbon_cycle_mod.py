"""
Module to handle carbon cycle from CO2 emissions to concentrations
"""

from functools import partial

import numpy as np
import pandas as pd

from .._utils import cut_and_check_pamset, update_pam_if_numeric
from ..pub_utils import _check_array_consistency
from .common_carbon_cycle_functions import (
    GE_COEFF,
    OCEAN_AREA,
    PPM_CO2_TO_PG_C,
    PPMKG_TO_UMOL_PER_VOL,
    PREINDUSTRIAL_CO2_CONC,
    calculate_airborne_fraction,
)
from .rfuns import rb_function, rb_function2, rs_function2, rs_function_array


def _process_flat_carbon_parameters(pamset):
    """
    Process flat carbon cycle parameters and convert to dictionary format.

    This function allows users to specify carbon cycle function parameters
    as individual floats (e.g., rb_coef0, rb_coef1, rb_tim0, rb_tim1)
    instead of requiring dictionary structures.

    Parameters
    ----------
    pamset : dict
        Parameter set that may contain flat carbon cycle parameters

    Returns
    -------
    dict
        Updated pamset with flat parameters converted to dictionary structures

    Examples
    --------
    Input: {"rb_coef0": 0.5, "rb_coef1": 0.25, "rb_tim0": 2.5, "rb_tim1": 10}
    Output: {"rb_function": {"coeffs": [0.5, 0.25], "timescales": [2.5, 10]}}
    """
    pamset = pamset.copy()

    # Process rb_function parameters
    rb_coefs = []
    rb_tims = []
    rb_keys_to_remove = []

    # Find all rb_coef* and rb_tim* parameters
    for key in pamset.keys():
        if key.startswith("rb_coef") and key[7:].isdigit():
            idx = int(key[7:])
            # Ensure list is long enough
            while len(rb_coefs) <= idx:
                rb_coefs.append(None)
            rb_coefs[idx] = pamset[key]
            rb_keys_to_remove.append(key)
        elif key.startswith("rb_tim") and key[6:].isdigit():
            idx = int(key[6:])
            # Ensure list is long enough
            while len(rb_tims) <= idx:
                rb_tims.append(None)
            rb_tims[idx] = pamset[key]
            rb_keys_to_remove.append(key)

    # Convert to rb_function dictionary if parameters found
    if rb_coefs or rb_tims:
        # Remove None values and check for consistency
        rb_coefs = [c for c in rb_coefs if c is not None]
        rb_tims = [t for t in rb_tims if t is not None]

        if len(rb_coefs) != len(rb_tims):
            raise ValueError(
                f"Number of rb_coef parameters ({len(rb_coefs)}) must match "
                f"number of rb_tim parameters ({len(rb_tims)})"
            )

        if len(rb_coefs) > 0:
            pamset["rb_function"] = {"coeffs": rb_coefs, "timescales": rb_tims}

        # Remove individual parameters
        for key in rb_keys_to_remove:
            pamset.pop(key)

    # Process rs_function parameters
    rs_coefs = []
    rs_tims = []
    rs_keys_to_remove = []

    # Find all rs_coef* and rs_tim* parameters
    for key in pamset.keys():
        if key.startswith("rs_coef") and key[7:].isdigit():
            idx = int(key[7:])
            # Ensure list is long enough
            while len(rs_coefs) <= idx:
                rs_coefs.append(None)
            rs_coefs[idx] = pamset[key]
            rs_keys_to_remove.append(key)
        elif key.startswith("rs_tim") and key[6:].isdigit():
            idx = int(key[6:])
            # Ensure list is long enough
            while len(rs_tims) <= idx:
                rs_tims.append(None)
            rs_tims[idx] = pamset[key]
            rs_keys_to_remove.append(key)

    # Convert to rs_function dictionary if parameters found
    if rs_coefs or rs_tims:
        # Remove None values and check for consistency
        rs_coefs = [c for c in rs_coefs if c is not None]
        rs_tims = [t for t in rs_tims if t is not None]

        # For rs_function, coeffs should have one more element than timescales
        if len(rs_coefs) != len(rs_tims) + 1:
            raise ValueError(
                f"For rs_function, number of rs_coef parameters ({len(rs_coefs)}) "
                f"must be one more than number of rs_tim parameters ({len(rs_tims)})"
            )

        if len(rs_coefs) > 0:
            pamset["rs_function"] = {"coeffs": rs_coefs, "timescales": rs_tims}

        # Remove individual parameters
        for key in rs_keys_to_remove:
            pamset.pop(key)

    return pamset


def sigmoid_gen(eval_point, sigmoid_center, sigmoid_width):
    """
    Generalised sigmoid function

    Parameters
    ----------
    eval_point : float
        The value at which to evaluate the sigmoid function
    sigmoid_center : float
        The center of the sigmoid function, where it is 1
    sigmoid_width : float
        The width of the sigmoid function

    Returns
    -------
        float
        sigmoid function value at eval_point
    """
    sigmoid = 1 / (1 + np.exp(-((eval_point - sigmoid_center) / sigmoid_width)))
    return sigmoid


def threshold_gen(eval_point, threshold_half, threshold_width):
    """
    Generalised threshold damping function

    Parameters
    ----------
    eval_point : float
        The value at which to evaluate the sigmoid function
    threshold_half : float
        The value for which the dampening of the threshold = 0.5
    threshold_width : float
        The width of the threshold

    Returns
    -------
        float
        threshold function value at eval_point
    """
    threshold = 1 - 1 / (
        1 + np.exp(-((eval_point - threshold_half) / threshold_width))
    )  # Threshold function
    return threshold


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


CARBON_CYCLE_MODEL_REQUIRED_PAMSET = {
    "beta_f": 0.287,
    "mixed_carbon": 75.0,
    "ml_w_sigmoid": 3.0,
    "ml_fracmax": 0.5,
    "ml_t_half": 0.5,
    "npp0": 60,
    "t_half": 0.5,
    "w_sigmoid": 7,
    "t_threshold": 4,
    "w_threshold": 7,
    "solubility_sens": 0.02,
    "solubility_limit": 0.5,
}


class CarbonCycleModel:
    """
    Class to handle carbon cycle calculations
    """

    def __init__(self, pamset_emiconc, pamset_carbon=None):
        """
        Initialise Carbon cycle model

        Parameters
        ----------
            pamset : dict
        """
        pamset = cut_and_check_pamset(
            {
                "idtm": 24,
                "nystart": 1750,
                "nyend": 2100,
            },
            pamset_emiconc,
        )
        if pamset_carbon is None:
            pamset_carbon = CARBON_CYCLE_MODEL_REQUIRED_PAMSET
        else:
            # Process flat carbon cycle parameters first
            pamset_carbon = _process_flat_carbon_parameters(pamset_carbon)
            pamset_carbon = cut_and_check_pamset(
                CARBON_CYCLE_MODEL_REQUIRED_PAMSET,
                pamset_carbon,
                used={"rs_function": "missing", "rb_function": "missing"},
            )

        pamset_carbon = take_out_missing(pamset_carbon.copy())
        self.pamset = {**pamset, **pamset_carbon}
        self.pamset["years_tot"] = pamset["nyend"] - pamset["nystart"] + 1
        self.reset_co2_hold(pamset_carbon)
        self.precalc_r_functions()

    def reset_co2_hold(self, pamset_carbon=None):
        """
        Reset values of CO2_hold for new run

        This method is mainly called to do a new run with the same cscm instance,
        in which case you need to reset hold values, and be able to update
        parameter values for the carbon cycle free parameter beta_f
        """
        self.co2_hold = {
            "yCO2": 0.0,
            "xCO2": PREINDUSTRIAL_CO2_CONC,
            "sCO2": np.zeros(self.pamset["idtm"] * self.pamset["years_tot"]),
            "emCO2_prev": 0.0,
            "dfnpp": np.zeros(self.pamset["idtm"] * self.pamset["years_tot"]),
            "ss1": 0.0,
            "sums": 0.0,
        }
        if pamset_carbon is not None:
            # Process flat carbon cycle parameters first
            pamset_carbon = _process_flat_carbon_parameters(pamset_carbon)
            # Check if we need to recompute r_functions
            needs_rfunction_recompute = (
                "rs_function" in pamset_carbon or "rb_function" in pamset_carbon
            )

            self.pamset = update_pam_if_numeric(
                self.pamset,
                pamset_new=pamset_carbon,
                can_change=CARBON_CYCLE_MODEL_REQUIRED_PAMSET.keys(),
            )
            # Update with non-numeric parameters (like function dictionaries)
            for key in ["rs_function", "rb_function"]:
                if key in pamset_carbon:
                    self.pamset[key] = pamset_carbon[key]
            self.fnpp_from_temp_vec = np.vectorize(self.fnpp_from_temp)

            # Recompute r_functions if we updated any function parameters
            if needs_rfunction_recompute:
                self.precalc_r_functions()

    def _set_co2_hold(
        self, xco2=PREINDUSTRIAL_CO2_CONC, yco2=0.0, emco2_prev=0.0, ss1=0.0, sums=0
    ):  # pylint: disable=too-many-positional-arguments, too-many-arguments
        """
        Reset the CO2 hold scalar values,

        Use this to rerun from same state as before
        in year not zero. Should only be used for back-calculations

        Parameters
        ----------
        xco2 : float
            CO2 concentration to set, default is PREINDUSTRIAL_CO2_CONC
            which is 278.0 andthe start value
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
            # self.pamset["rs_function"] = rs_function_array
            self.r_functions[0, :] = rs_function_array(
                np.arange(self.pamset["idtm"] * self.pamset["years_tot"]),
                self.pamset["idtm"],
            )
        else:
            coeffs, timscales = _check_array_consistency(
                self.pamset["rs_function"]["coeffs"],
                self.pamset["rs_function"]["timescales"],
                for_rs=True,
            )
            self.r_functions[0, :] = rs_function2(
                np.arange(self.pamset["idtm"] * self.pamset["years_tot"]),
                rs_coef=coeffs,
                rs_tim=timscales,
                idtm=self.pamset["idtm"],
            )

        if "rb_function" not in self.pamset:
            self.r_functions[1, :] = rb_function(
                np.arange(self.pamset["idtm"] * self.pamset["years_tot"]),
                self.pamset["idtm"],
            )
        else:
            coeffs, timscales = _check_array_consistency(
                self.pamset["rb_function"]["coeffs"],
                self.pamset["rb_function"]["timescales"],
            )
            self.r_functions[1, :] = rb_function2(
                np.arange(self.pamset["idtm"] * self.pamset["years_tot"]),
                rb_coef=coeffs,
                rb_tim=timscales,
                idtm=self.pamset["idtm"],
            )

    def fnpp_from_temp(self, dtemp=0):
        """
        Temperature dependence function for fnpp

        Parameters
        ----------
        dtemp : float
            Degrees of temperature since start of run

        Returns
        -------
        float
            fnpp at given temperature for assumed linear
            relationship
        """
        sigmoid = partial(
            sigmoid_gen,
            sigmoid_center=self.pamset["t_half"],
            sigmoid_width=self.pamset["w_sigmoid"],
        )
        threshold = partial(
            threshold_gen,
            threshold_half=self.pamset["t_threshold"],
            threshold_width=self.pamset["w_threshold"],
        )
        # Threshold function
        baseline = self.pamset["npp0"] / (sigmoid(0) * threshold(0))

        return baseline * sigmoid(dtemp) * threshold(dtemp)

    def solubility_temp_feedback(self, dtemp=0.0):
        """
        Exponential scaling of CO2 solubility with temperature, with upper limit.

        The parameters used are solubility_sens, the fractional decrease in solubility
        per degree C (default 0.02) and solubility_limit, the maximum allowed gain
        (e.g. 0.5 for 50% increase)

        Parameters
        ----------
        dtemp : float
            Degrees of temperature since start of run


        Returns
        -------
        float
            Solubility scaling factor (max 1 + solubility_limit)
        """
        scale = np.exp(-self.pamset["solubility_sens"] * dtemp)
        scale = np.clip(scale, 0, 1 + self.pamset["solubility_limit"])
        # Limit the gain to 1 + solubility_limit (e.g. 1.5)
        return scale

    def mixed_layer_temp_feedback(self, dtemp=0):
        """
        Calculate mixed layer depth change with temperature using new parameters.

        Parameters
        ----------
        dtemp : float
            Degrees of temperature since start of run

        Returns
        -------
        float
            Mixed layer depth change factor
        """
        sigmoid = partial(
            sigmoid_gen,
            sigmoid_center=self.pamset["ml_t_half"],
            sigmoid_width=self.pamset["ml_w_sigmoid"],
        )
        mixed_layer_temp = 1 - self.pamset["ml_fracmax"] * sigmoid(dtemp)
        baseline = 1 - self.pamset["ml_fracmax"] * sigmoid(0)
        return self.pamset["mixed_carbon"] * mixed_layer_temp / baseline

    def _calculate_partial_pressure_mixed_layer(self, it, dtemp=0):
        """
        Calculate ocean mixed layer partial pressure

        Parameters
        ----------
        it : int
            Subyearly timestep at which to calculate
        dtemp : float
            Temperature change at timestep

        Returns
        -------
            float
            Mixed ocean layer partial pressure
        """
        dt = 1 / self.pamset["idtm"]
        if it > 0:
            # Pulse response integrate carbon content in the mixed layer
            # Carbon decays into deep layer according to pulse response function
            sumz = np.dot(
                self.co2_hold["sCO2"][: it - 1], np.flip(self.r_functions[0, 1:it])
            )
        else:
            sumz = 0.0

        # Inorganic carbon content in the mixed layer:
        z_co2 = (
            PPMKG_TO_UMOL_PER_VOL
            * GE_COEFF
            * dt
            / self.mixed_layer_temp_feedback(dtemp)
            * (sumz + 0.5 * self.co2_hold["sCO2"][it])
        )
        # Partial pressure in ocean mixed layer,
        # in principle this only holds in 17.7-18.2 temperature range
        # but up to 1320 ppm
        # and this particular formulation is the solution at T = 18.2
        # The original description paper also has a different formulation
        # with a wider valid temperature range, but only up to 200 ppm
        # which is not used
        # This might be a natural place to look for/substitute with a
        # different / more general / temperature dependent formulation
        # which would be in line with the model philosophy and structure
        return self.solubility_temp_feedback(dtemp) * (
            1.3021 * z_co2
            + 3.7929e-3 * (z_co2**2)
            + 9.1193e-6 * (z_co2**3)
            + 1.488e-8 * (z_co2**4)
            + 1.2425e-10 * (z_co2**5)
        )

    def co2em2conc(
        self, yr, em_co2_common, dtemp=0.0
    ):  # pylint: disable=too-many-locals
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
        dtemp : float
            temperature change from start of run at previous timestep

        Returns
        -------
        float
             CO2 concetrations for year in question
        """
        # TIMESTEP (YR)
        dt = 1.0 / self.pamset["idtm"]

        cc1 = dt * OCEAN_AREA * GE_COEFF / (1 + dt * OCEAN_AREA * GE_COEFF / 2.0)
        yr_ix = yr - self.pamset["nystart"]
        fnpp = self.fnpp_from_temp(dtemp=dtemp)
        # Monthloop:
        for i in range(self.pamset["idtm"]):
            it = yr_ix * self.pamset["idtm"] + i
            sumf = 0.0

            # Net emissions, including biogenic fertilization effects
            if it > 0:
                # Net primary production in timestep
                self.co2_hold["dfnpp"][it] = (
                    fnpp
                    * self.pamset["beta_f"]
                    * np.log(self.co2_hold["xCO2"] / PREINDUSTRIAL_CO2_CONC)
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
                # Trapesoidal part of air-sea flux integral for first timestep
                self.co2_hold["ss1"] = 0.5 * em_co2 / (OCEAN_AREA * GE_COEFF)
                ss2 = self.co2_hold["ss1"]
                self.co2_hold["sums"] = 0.0
            else:
                # Trapesoidal part of air-sea flux integral for last timestep
                ss2 = 0.5 * em_co2 / (OCEAN_AREA * GE_COEFF) - self.co2_hold["yCO2"] / (
                    dt * OCEAN_AREA * GE_COEFF
                )
                # Integral inner part for air-sea flux
                self.co2_hold["sums"] = (
                    self.co2_hold["sums"]
                    + self.co2_hold["emCO2_prev"] / (OCEAN_AREA * GE_COEFF)
                    - self.co2_hold["sCO2"][it - 1]
                )
            # Air-sea carbon flux adding the integral parts for the trapezoidal rule
            self.co2_hold["sCO2"][it] = cc1 * (
                self.co2_hold["sums"] + self.co2_hold["ss1"] + ss2
            )
            self.co2_hold["emCO2_prev"] = em_co2
            # Partial pressure in ocean mixed layer,
            # in principle this only holds in 17.7-18.2 temperature range
            # but up to 1320 ppm
            # and this particular formulation is the solution at T = 18.2
            # The original description paper also has a different formulation
            # with a wider valid temperature range, but only up to 200 ppm
            # which is not used
            # This might be a natural place to look for/substitute with a
            # different / more general / temperature dependent formulation
            # which would be in line with the model philosophy and structure
            yco2_prev = self.co2_hold["yCO2"]
            self.co2_hold["yCO2"] = self._calculate_partial_pressure_mixed_layer(it)
            # Partial pressure in the atmosphere, this comes from
            # solving the transfer equation between atmosphere and
            # ocean  to get the resulting atmosphere partial pressure
            # We use the midpoint partial pressure here as this ensures
            # best (but not perfect) closure of the carbon cycle
            self.co2_hold["xCO2"] = (
                self.co2_hold["sCO2"][it]
                + 0.5 * self.co2_hold["yCO2"]
                + 0.5 * yco2_prev
                + PREINDUSTRIAL_CO2_CONC
            )
            # print("it: %d, emCO2: %e, sCO2: %e, zCO2: %e, yCO2: %e, xCO2: %e, ss1: %e, ss2: %e, dnfpp:%e"%(it, em_co2, self.co2_hold["sCO2"][it], z_co2, self.co2_hold["yCO2"], self.co2_hold["xCO2"], self.co2_hold["ss1"], ss2, self.co2_hold["dfnpp"][it]))
        return self.co2_hold["xCO2"]

    def _get_ffer_timeseries(
        self, conc_run=False, co2_conc_series=None, dtemp_series=None
    ):
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
            if dtemp_series is None:
                dtemp_series = np.zeros(len(co2_conc_series))
            timesteps = len(co2_conc_series) * self.pamset["idtm"]
            fnpp = np.repeat(
                self.fnpp_from_temp_vec(
                    dtemp=dtemp_series,
                ),
                self.pamset["idtm"],
            )
            dfnpp = (
                np.repeat(
                    [
                        self.pamset["beta_f"]
                        * np.log(co2_conc / PREINDUSTRIAL_CO2_CONC)
                        for co2_conc in co2_conc_series
                    ],
                    self.pamset["idtm"],
                )
                * fnpp
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

    def get_biosphere_carbon_flux(
        self, conc_run=False, co2_conc_series=None, dtemp_series=None
    ):
        """
        Get the carbon flux to the biosphere as timeseries over years

        For emissions runs this yields appropriate amounts only for years for which the
        model has been run and zeros otherwise. If a concentrations series is used
        this yields a back calculated estimate of the carbon flux the model
        estimates given these atmospheric concentrations.

        Parameters
        ----------
        conc_run : bool
            Whether calculation is being done for a conc_run
            (i.e. a run that has CO2_hold values precalculated, or whether back
            calculations are needed or co2_conc_series need to be used for calculations)
        co2_conc_series : np.ndarray
            Time series of concentrations to calculate ffer from
        dtemp_series : np.ndarray
            Timeseries of temperature change, should be included for calculation of
            outputs if temperature feedbacks are on

        Returns
        -------
        np.ndarray
            Timeseries of the yearly carbon flux to the biosphere
        """
        ffer = self._get_ffer_timeseries(conc_run, co2_conc_series, dtemp_series)
        biosphere_carbon_flux = (
            np.array(
                [
                    np.sum(
                        ffer[
                            self.pamset["idtm"]
                            * yrix : self.pamset["idtm"]
                            * (yrix + 1)
                        ]
                    )
                    for yrix in range(self.pamset["years_tot"])
                ]
            )
            / self.pamset["idtm"]
        )

        return biosphere_carbon_flux

    def get_ocean_carbon_flux(
        self, conc_run=False, co2_conc_series=None, dtemp_series=None
    ):
        """
        Get yearly timeseries of ocean carbon flux

        If this is called from a concentrations either it is
        assumed that emissions are already back calculated
        Otherwhise you can send a concentration series to do
        the back calculations with to set all the values in
        the sCO2 array

        Returns
        -------
        np.ndarray
            Timeseries of the added carbon content to the yearly
            ocean carbon flux (Pg / C /yr)
        """
        # TODO: Don't redo this if back-calculation is already done...
        if conc_run and co2_conc_series is not None:
            self.back_calculate_emissions(co2_conc_series, dtemp_series=dtemp_series)
        ocean_carbon_flux = (
            np.array(
                [
                    np.sum(
                        self.co2_hold["sCO2"][
                            self.pamset["idtm"]
                            * yrix : self.pamset["idtm"]
                            * (yrix + 1)
                        ]
                    )
                    for yrix in range(self.pamset["years_tot"])
                ]
            )
            * PPM_CO2_TO_PG_C
            / self.pamset["idtm"]
            * GE_COEFF
            * OCEAN_AREA
        )
        return ocean_carbon_flux

    def back_calculate_emissions(self, co2_conc_series, dtemp_series=None):
        """
        Back calculate emissions from conc_run

        co2_conc_series is assumed to be the series of concentrations
        from the year 0

        Parameters
        ----------
        co2_conc_series : np.ndarray
            Timeseries of co2 concentrations for which to back
            calculate emissions
        dtemp_series : np.ndarray
            Timeseries of temperature change for which to back calculate emissions
            It should be the same length as the concentration timeseries
            If no value is sent, a timeseries of zeros will be used

        Returns
        -------
            np.ndarray
            Timeseries of estimated emissions to match the concentration and
            temperature timeseries sent
        """
        prev_co2_conc = PREINDUSTRIAL_CO2_CONC
        em_series = np.zeros(len(co2_conc_series))
        if dtemp_series is None:
            dtemp_series = np.zeros(len(co2_conc_series))

        ffer = self._get_ffer_timeseries(
            conc_run=True,
            co2_conc_series=co2_conc_series,
            dtemp_series=dtemp_series,
        )
        for i, co2_conc in enumerate(co2_conc_series):
            ffer_here = ffer[i * self.pamset["idtm"]]
            em_series[i] = self._guess_emissions_iteration(
                co2_conc,
                prev_co2_conc,
                yrix=i,
                ffer=ffer_here,
                dtemp=dtemp_series[i],
            )
            prev_co2_conc = co2_conc
        return em_series

    def _guess_emissions_iteration(
        self,
        co2_conc_now,
        co2_conc_zero,
        dtemp=0,
        yrix=0,
        rtol=1e-7,
        maxit=100,
        ffer=None,
    ):  # pylint: disable=too-many-arguments, too-many-positional-arguments
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
            ffer = self._get_ffer_timeseries(
                [co2_conc_zero, co2_conc_now], dtemp_series=[0, dtemp]
            )[yrix * self.pamset["idtm"]]
        co2_change = co2_conc_now - co2_conc_zero
        min_guess = np.min(
            (
                co2_change * PPM_CO2_TO_PG_C + 8 * ffer,
                co2_change * PPM_CO2_TO_PG_C - 4 * ffer,
            )
        )  # self.simplified_em_backward(co2_conc_zero - 4*co2_change , co2_conc_zero)
        max_guess = np.max(
            (
                co2_change * PPM_CO2_TO_PG_C + 8 * ffer,
                co2_change * PPM_CO2_TO_PG_C - 4 * ffer,
            )
        )  # self.simplified_em_backward(co2_conc_zero + 4* co2_change, co2_conc_zero)
        if max_guess - min_guess < 1:
            max_guess = max_guess + 1
            min_guess = min_guess - 1
        guess = np.mean(
            (min_guess, max_guess)
        )  # self.simplified_em_backward(co2_conc_now, co2_conc_zero)
        hold_dict = self._get_co2_hold_values()
        estimated_conc = self.co2em2conc(
            self.pamset["nystart"] + yrix, guess, dtemp=dtemp
        )
        iteration = 0
        if yrix % 50 == 0:
            print(f"yr: {yrix} has minguess: {min_guess}, maxguess: {max_guess}")
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
                self.pamset["nystart"] + yrix, guess, dtemp=dtemp
            )
            iteration = iteration + 1
        if yrix % 50 == 0:
            print(
                f"End guess: {guess} {co2_conc_now} and {co2_conc_zero}  and estimated conc {estimated_conc}"
            )
        return guess

    def get_carbon_cycle_output(
        self, years, conc_run=False, conc_series=None, dtemp_series=None
    ):
        """
        Make and return a dataframe with carbon cycle data

        Parameters
        ----------
        years : np.array
            Array of years to use as index
        conc_run : bool
            Whether this is from a concentrations driven run or emissions driven
        conc_series : np.array
            Numpy array of concentrations, must be included for a concentrations
            driven run to back calculate emissions
        dtemp_series : np.ndarray
            Timeseries of temperature change, should be included for calculation of
            outputs if temperature feedbacks are on

        Returns
        -------
        Pandas.DataFrame
            Including carbon cycle variables
        """
        if conc_run and conc_series is None:
            return None
        if dtemp_series is None:
            dtemp_series = np.zeros(self.pamset["years_tot"])
        if conc_run:
            em_series = self.back_calculate_emissions(conc_series, dtemp_series)
        df_carbon = pd.DataFrame(
            data={
                "Biosphere carbon flux": self.get_biosphere_carbon_flux(
                    conc_run=conc_run,
                    co2_conc_series=conc_series,
                    dtemp_series=dtemp_series,
                ),
                "Ocean carbon flux": self.get_ocean_carbon_flux(),
                "Mixed layer depth": self.mixed_layer_temp_feedback(dtemp_series),
                "CO2 solubility": self.solubility_temp_feedback(dtemp_series),
                "Net Primary Production": self.fnpp_from_temp_vec(dtemp_series),
            },
            index=years,
        )

        if not conc_run:
            return df_carbon
        df_carbon["Airborne fraction CO2"] = calculate_airborne_fraction(
            em_series, conc_series  # pylint: disable=possibly-used-before-assignment
        )
        df_carbon["Emissions"] = em_series
        return df_carbon
