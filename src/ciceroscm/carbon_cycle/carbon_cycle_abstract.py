"""
Abstract carbon cycle model
"""

from abc import ABC, abstractmethod

import numpy as np

from .._utils import update_pam_if_numeric
from .common_carbon_cycle_functions import (
    PPM_CO2_TO_PG_C,
    PREINDUSTRIAL_CO2_CONC,
    carbon_cycle_init_pamsets,
)


class AbstractCarbonCycleModel(ABC):
    """
    Abstract class to define carbon cycle methodology
    """

    @property
    @abstractmethod
    def carbon_cycle_model_required_pamset(self):
        """
        Overwrite this to set a required carbon cycle
        model parameterset dictionary. This should be
        a class variable
        """

    @classmethod
    def get_carbon_cycle_required_pamset(cls):
        """
        Get the class variable carbon_cycle_model_required_pamset
        which should be a class variable. This setup is to
        substitute for the lack of abstract class variables in python
        """
        return cls.carbon_cycle_model_required_pamset

    def __init__(self, pamset_emiconc, pamset_carbon):
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
            Optional carbon specific parameter set allowed options should be defined
            in the inherited class
        """
        pamset, pamset_carbon = carbon_cycle_init_pamsets(
            pamset_emiconc, pamset_carbon, self.get_carbon_cycle_required_pamset()
        )
        self.pamset = {**pamset, **pamset_carbon}
        self.pamset["years_tot"] = pamset["nyend"] - pamset["nystart"] + 1
        self.reset_co2_hold(pamset_carbon)

    def reset_co2_hold(self, pamset_carbon=None):
        """
        Reset values of CO2_hold for new run

        This method is mainly called to do a new run with the same cscm instance,
        in which case you might need to reset hold values
        (this needs to be implemented concretely in subclass),
        and be able to update parameter values for the carbon cycle free parameters

        Parameters
        ----------
        pamset_carbon : dict
            Optional dictionary of new values for a subset of the free
            carbon cycle parameters
        """
        if pamset_carbon is not None:
            self.pamset = update_pam_if_numeric(
                self.pamset,
                pamset_new=pamset_carbon,
                can_change=self.get_carbon_cycle_required_pamset().keys(),
            )

    @abstractmethod
    def co2em2conc(self, yr, em_co2_common, feedback_dict=None):
        """
        Calculate co2 concentrations from emissions

        Method to calculate co2 concentrations from emissions
        Needs to be implemented in subclass

        Parameters
        ----------
        yr : int
          Year for which to calculate
        em_co2_common : float
             Sum of CO2 emissions from fossil fuels, land use change and natural emissions
             for the year in question
        feedback_dict : dict
            Dictionary containing feedback variables and their values

        Returns
        -------
        float
             CO2 concetrations for year in question
        """

    def _get_co2_hold_values(self):
        """
        Get a dictionary of state variables for the carbon cycle that allows you
        to reset your model state, typically pool sizes, chemical composition
        current year etc...

        Make sure to overwrite this is if your model actually has state dependence
        """
        return {}

    def _set_co2_hold(self, hold_dict=None):
        """
        Set state variables for the carbon cycle from hold_dict dictionary
        thereby resetting your model state, typically pool sizes, chemical composition
        current year etc...

        Make sure to overwrite this is if your model actually has state dependence
        """

    def get_feedback_list(self):
        """
        Get a list of feedback variables for the carbon cycle model

        Overwrite this if your carbon cycle model has needs more or other feedback variables
        If your carbon cycle doesn't need any feedback variables, you should return None

        Returns
        -------
        list or None
            List of feedback variable names or None if no feedback variables are needed
        """
        return ["dtemp"]

    def _guess_emissions_iteration(
        self,
        co2_conc_now,
        initial_max_min_guess,
        feedback_dict=None,
        yrix=0,
        rtol=1e-7,
        maxit=100,
    ):  # pylint: disable=too-many-arguments, too-many-positional-arguments
        """
        Iterate to get right emissions for a single year

        Make sure the yrix is where you are at

        Parameters
        ----------
        co2_conc_now : float
            Value of CO2 concentration in timestep resulting after
            the emissions you would like to find are applied
        initial_max_min_guess : tuple
            Initial maximum and minimum guess for emissions
            to start the bisection iteration
        feedback_dict : dict
            Dictionary containing feedback variables and their values
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

        Returns
        -------
        float
            Back calculated emissions that yields a concentrations change
            from the conc_co2_zero to conc_co2_now which

        """
        max_guess = initial_max_min_guess[0]
        min_guess = initial_max_min_guess[1]
        guess = np.mean((min_guess, max_guess))
        hold_dict = self._get_co2_hold_values()
        estimated_conc = self.co2em2conc(
            self.pamset["nystart"] + yrix, guess, feedback_dict=feedback_dict
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
            self._set_co2_hold(hold_dict)
            estimated_conc = self.co2em2conc(
                self.pamset["nystart"] + yrix, guess, feedback_dict=feedback_dict
            )
            iteration = iteration + 1
        return guess

    def get_initial_max_min_guess(
        self, co2_conc_now, co2_conc_zero, yrix=0, feedback_dict=None, em_width=None
    ):  # pylint: disable=too-many-positional-arguments, too-many-arguments, unused-argument
        """
        Calculate initial max and min guess i.e.
        maximum span of emissions that could have lead to the
        change from co2_conc_zero to co2_conc_now
        I.e. starting guess span for bisection
        """
        if feedback_dict is None:
            feedback_dict_series = {
                key: np.zeros(2) for key in self.get_feedback_list()
            }
        else:
            feedback_dict_series = {
                key: np.array([feedback_dict.get(key, 0), feedback_dict.get(key, 0)])
                for key in feedback_dict
            }
        co2_change = co2_conc_now - co2_conc_zero
        if em_width is None:
            em_width = self.make_guess_emsize_estimates(
                [co2_conc_zero, co2_conc_now], feedback_dict_series=feedback_dict_series
            )[1]

        min_guess = np.min(
            (
                co2_change * PPM_CO2_TO_PG_C + 8 * em_width,
                co2_change * PPM_CO2_TO_PG_C - 4 * em_width,
            )
        )  # self.simplified_em_backward(co2_conc_zero - 4*co2_change , co2_conc_zero)
        max_guess = np.max(
            (
                co2_change * PPM_CO2_TO_PG_C + 8 * em_width,
                co2_change * PPM_CO2_TO_PG_C - 4 * em_width,
            )
        )  # self.simplified_em_backward(co2_conc_zero + 4* co2_change, co2_conc_zero)
        if max_guess - min_guess < 1:
            max_guess = max_guess + 1
            min_guess = min_guess - 1

        return max_guess, min_guess

    def make_guess_emsize_estimates(
        self, co2_conc_series, feedback_dict_series=None
    ):  # pylint: disable=unused-argument
        """
        Make an estimate of emission size changes from a concentrations
        change and a temperature size. This should just be used to scale
        the possible emissions size for the bisection in backcalculation

        This method should mostly be overwritten by concrete carbon cycle
        as information about pool states or impulse response historical /
        hysteresis effects might need to be taken into account

        Parameters
        ----------
        co2_conc_series : np.ndarray
            Timeseries of co2 concentrations for which to back
            calculate emissions
        feedback_dict_series : np.ndarray
            Timeseries of feedback variables for which to back calculate emissions
            It should be the same length as the concentration timeseries
            If no value is sent, a timeseries of zeros will be used
            In this abstract class we only do not use this variable
            but it might be useful in concrete implementations

        Returns
        -------
            np.ndarray
            Timeseries of estimated emissions very rough size of possible
            emission changes that can be used
        """
        prev_co2_series = np.zeros_like(co2_conc_series)
        prev_co2_series[0] = PREINDUSTRIAL_CO2_CONC
        prev_co2_series[1:] = co2_conc_series[:-1]
        return (co2_conc_series - prev_co2_series) * PPM_CO2_TO_PG_C

    def back_calculate_emissions(self, co2_conc_series, feedback_dict_series=None):
        """
        Back calculate emissions from conc_run

        co2_conc_series is assumed to be the series of concentrations
        from the year 0

        Parameters
        ----------
        co2_conc_series : np.ndarray
            Timeseries of co2 concentrations for which to back
            calculate emissions
        feedback_dict_series : dict
            Timeseries of feedback variables for which to back calculate emissions
            For each variable key, the value should be an np.ndarray timeseries of
            the same length as the concentration timeseries
            If no value is sent, timeseries of zeros will be used

        Returns
        -------
            np.ndarray
            Timeseries of estimated emissions to match the concentration and
            temperature timeseries sent
        """
        prev_co2_conc = PREINDUSTRIAL_CO2_CONC
        em_series = np.zeros(len(co2_conc_series))
        if feedback_dict_series is None:
            feedback_dict_series = {
                key: np.zeros(len(co2_conc_series)) for key in self.get_feedback_list()
            }
        emissions_width = self.make_guess_emsize_estimates(
            co2_conc_series=co2_conc_series,
            feedback_dict_series=feedback_dict_series,
        )
        for i, co2_conc in enumerate(co2_conc_series):
            em_width_here = emissions_width[i * self.pamset["idtm"]]
            em_series[i] = self._guess_emissions_iteration(
                co2_conc,
                self.get_initial_max_min_guess(
                    co2_conc, prev_co2_conc, em_width=em_width_here
                ),
                yrix=i,
                feedback_dict={
                    key: feedback_dict_series[key][i] for key in feedback_dict_series
                },
            )
            prev_co2_conc = co2_conc
        return em_series

    @abstractmethod
    def get_carbon_cycle_output(
        self, years, conc_run=False, conc_series=None, feedback_dict_series=None
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
        feedback_dict_series : dict
            Timeseries of feedback variables for which to back calculate emissions
            For each variable key, the value should be an np.ndarray timeseries of
            the same length as the concentration timeseries
            If no value is sent, timeseries of zeros will be used

        Returns
        -------
        Pandas.DataFrame
            Including carbon cycle variables
        """
