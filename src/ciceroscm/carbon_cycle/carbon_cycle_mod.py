"""
Module to handle carbon cycle from CO2 emissions to concentrations
"""

from functools import partial

import numpy as np
import pandas as pd

from .._utils import update_pam_if_numeric
from ..pub_utils import _check_array_consistency
from .carbon_cycle_abstract import AbstractCarbonCycleModel
from .common_carbon_cycle_functions import (
    GE_COEFF,
    OCEAN_AREA,
    PPM_CO2_TO_PG_C,
    PPMKG_TO_UMOL_PER_VOL,
    PREINDUSTRIAL_CO2_CONC,
    calculate_airborne_fraction,
    carbon_cycle_init_pamsets,
)
from .rfuns import (
    _process_flat_carbon_parameters,
    rb_function,
    rb_function2,
    rs_function2,
    rs_function_array,
)


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


CARBON_CYCLE_MODEL_REQUIRED_PAMSET = {
    "beta_f": 0.287,
    "mixed_carbon": 75.0,
    "ml_w_sigmoid": 3.0,
    "ml_fracmax": 0.5,
    "ml_t_half": 0.5,
    "t_half": 0.5,
    "w_sigmoid": 7,
    "t_threshold": 4,
    "w_threshold": 7,
    "solubility_sens": 0.02,
    "solubility_limit": 0.5,
}


class CarbonCycleModel(AbstractCarbonCycleModel):
    """
    Class to handle carbon cycle calculations
    """

    carbon_cycle_model_required_pamset = CARBON_CYCLE_MODEL_REQUIRED_PAMSET

    def __init__(
        self, pamset_emiconc, pamset_carbon=None
    ):  # pylint: disable=super-init-not-called
        """
        Initialise Carbon cycle model

        Parameters
        ----------
        pamset_emiconc : dict
            Parameter set from the concentrations emission handler, containing:
            - idtm: Number of subyearly timesteps (e.g., 24 for monthly steps).
            - nystart: Start year of the simulation.
            - nyend: End year of the simulation.
        pamset_carbon : dict
            Optional carbon specific parameter set allowed options
            - beta_f: CO2 fertilization factor (affects land carbon uptake).
            - mixed_carbon: depth of ocean mixed layer as seen by the carbon cycle (m)
            - ml_w_sigmoid: Mixed layer dampening sigmoid width (K).
            - ml_fracmax: Mixed layer max fractional dampening (0-1 unitless)
            - ml_t_half: Temperature change of mixed layer dampening sigmoid center (K)
            - t_half: Centerpoint of land uptake sigmoid (K),
            - w_sigmoid: Width of land uptake sigmoid (K),
            - t_threshold: Threshold temperature for land uptake decay (K),
            - w_threshold: Width of land uptake decay (K),
            - solubility_sen: Sensitivity of solubility of CO2 in the ocean to temperature.
            - solubility_limit: Limit to temperature solubility gain with temperature.
            - rs_function: User defined mixed layer to deep ocean impulse response function
            - rb_function: User defined land carbon impulse response decay function
        """
        if pamset_carbon is not None:
            pamset_carbon = _process_flat_carbon_parameters(pamset_carbon)
        pamset, pamset_carbon = carbon_cycle_init_pamsets(
            pamset_emiconc,
            pamset_carbon,
            self.get_carbon_cycle_required_pamset(),
            used={"rs_function": "missing", "rb_function": "missing"},
        )
        self.pamset = {**pamset, **pamset_carbon}
        self.pamset["years_tot"] = pamset["nyend"] - pamset["nystart"] + 1
        self.reset_co2_hold(pamset_carbon)
        self.precalc_r_functions()

    def reset_co2_hold(self, pamset_carbon=None):
        """
        Reset values of CO2_hold for new run

        This method is mainly called to do a new run with the same cscm instance,
        in which case you need to reset hold values, and be able to update
        parameter values for the carbon cycle free parameters

        Parameters
        ----------
        pamset_carbon : dict
            Optional dictionary of new values for a subset of the free
            carbon cycle parameters
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
        self, hold_dict=None
    ):  # pylint: disable=too-many-positional-arguments, too-many-arguments
        """
        Reset the CO2 hold scalar values,

        Use this to rerun from same state as before
        in year not zero. Should only be used for back-calculations

        Parameters
        ----------
        hold_dict : dict
            Dictionary containing hold values to set, the following keys are used.
            Defaults equal values for starting a run from pre-industrial conditions
            xco2 : float
                CO2 concentration to set, default is PREINDUSTRIAL_CO2_CONC
                which is 278.0 and the start value
            yco2 : float
                yco2 value, default is 0.0 which is the start value
            emco2_prev : float
                emissions in previous timestep, default is 0.0
            ss1 : float
                ss1 value, default is 0.0
            sums : float
                sums of ocean uptake inorganic carbon, default is 0.0
        """
        if hold_dict is None:
            hold_dict = {}
        self.co2_hold["yCO2"] = hold_dict.get("yco2", 0.0)
        self.co2_hold["xCO2"] = hold_dict.get("xco2", PREINDUSTRIAL_CO2_CONC)
        self.co2_hold["emCO2_prev"] = hold_dict.get("emco2_prev", 0.0)
        self.co2_hold["ss1"] = hold_dict.get("ss1", 0.0)
        self.co2_hold["sums"] = hold_dict.get("sums", 0.0)

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
        # Threshold function, 60. is initial npp-value
        baseline = 60.0 / (sigmoid(0) * threshold(0))

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
        self, yr, em_co2_common, feedback_dict=None
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
        feedback_dict : dict
            Dictionary containing feedback information, in this case key "dtemp"
            with value temperature change from pre-industrial in degrees K

        Returns
        -------
        float
             CO2 concetrations for year in question
        """
        if feedback_dict is None:
            dtemp = 0.0
        else:
            dtemp = feedback_dict.get("dtemp", 0.0)
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
            self.back_calculate_emissions(
                co2_conc_series, feedback_dict_series={"dtemp": dtemp_series}
            )
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

    def make_guess_emsize_estimates(self, co2_conc_series, feedback_dict_series=None):
        """
        Make an estimate of emission size changes from a concentrations
        change and a temperature size. This should just be used to scale
        the possible emissions size for the bisection in backcalculation

        Here we make this estimate by using the ffer timeseries as an
        assumptions.

        Parameters
        ----------
        co2_conc_series : np.ndarray
            Timeseries of co2 concentrations for which to back
            calculate emissions
        feedback_dict_series : dict
            Dictionary containing feedback variables and their values
            If no value is sent, a dictionary of zeros will be used
            For this carbon cycle model the only relevant variable is
            'dtemp' which is the temperature change at each timestep

        Returns
        -------
            np.ndarray
            Timeseries of estimated emissions very rough size of possible
            emission changes that can be used
        """
        if feedback_dict_series is None:
            return self._get_ffer_timeseries(
                conc_run=True,
                co2_conc_series=co2_conc_series,
                dtemp_series=np.zeros(len(co2_conc_series)),
            )
        return self._get_ffer_timeseries(
            conc_run=True,
            co2_conc_series=co2_conc_series,
            dtemp_series=feedback_dict_series.get(
                "dtemp", np.zeros(len(co2_conc_series))
            ),
        )

    def get_carbon_cycle_output(
        self, years, conc_run=False, conc_series=None, feedback_dict_series=None
    ):  # pylint: disable=unused-argument
        """
        Make and return a dataframe with carbon cycle data

        Parameters
        ----------
        years : np.array
            Array of years to use as index
        conc_run : bool
            Whether this is from a concentrations driven run or emissions driven
        conc_series : np.array
            Numpy array of concentrations of same length as years, m
            ust be included for a concentrations driven run to back calculate emissions
        feedback_dict_series : dict
            Dictionary containing feedback variables and their values
            For each variable key this should be a numpy array of same length as years
            If no value is sent, a dictionary of zeros will be used.
            This will yield wrong values if temperature feedbacks are on.
            For this carbon cycle model the only relevant variable is
            'dtemp' which is the temperature change at each timestep

        Returns
        -------
        Pandas.DataFrame
            Including carbon cycle variables
        """
        if conc_run and conc_series is None:
            return None
        if feedback_dict_series is None:
            feedback_dict_series = {"dtemp": np.zeros(len(years))}
        if conc_run:
            em_series = self.back_calculate_emissions(
                conc_series, feedback_dict_series=feedback_dict_series
            )
        df_carbon = pd.DataFrame(
            data={
                "Biosphere carbon flux": self.get_biosphere_carbon_flux(
                    conc_run=conc_run,
                    co2_conc_series=conc_series,
                    dtemp_series=feedback_dict_series.get("dtemp", None),
                ),
                "Ocean carbon flux": self.get_ocean_carbon_flux(
                    conc_run=conc_run,
                    co2_conc_series=conc_series,
                    dtemp_series=feedback_dict_series.get("dtemp", None),
                ),
                "Mixed layer depth": self.mixed_layer_temp_feedback(
                    feedback_dict_series.get("dtemp", None)
                ),
                "CO2 solubility": self.solubility_temp_feedback(
                    feedback_dict_series.get("dtemp", None)
                ),
                "Net Primary Production": self.fnpp_from_temp_vec(
                    feedback_dict_series.get("dtemp", None)
                ),
            },
            index=years,
        )

        if not conc_run:
            return df_carbon
        df_carbon["Airborne fraction CO2"] = calculate_airborne_fraction(
            em_series,  # pylint: disable=possibly-used-before-assignment
            conc_series,
        )
        df_carbon["Emissions"] = em_series
        return df_carbon
