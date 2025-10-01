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

    # TODO: Generalise so more than just temperature can be sent from thermal
    @abstractmethod
    def co2em2conc(self, yr, em_co2_common, dtemp=0.0):
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
        dtemp : float
            temperature change from start of run at previous timestep

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

    def _set_co2_hold_values(self, hold_dict):
        """
        Set state variables for the carbon cycle from hold_dict dictionary
        thereby resetting your model state, typically pool sizes, chemical composition
        current year etc...

        Make sure to overwrite this is if your model actually has state dependence
        """

    def _guess_emissions_iteration(
        self,
        co2_conc_now,
        co2_conc_zero,
        initial_max_min_guess,
        dtemp=0,
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
        max_guess = initial_max_min_guess[0]
        min_guess = initial_max_min_guess[1]
        guess = np.mean((min_guess, max_guess))
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

    @abstractmethod
    def get_initial_max_min_guess(
        self, co2_conc_now, co2_conc_zero, yrix=0, dtemp=0, ffer=None
    ):
        """
        Calculate initial max and min guess i.e.
        maximum span of emissions that could have lead to the
        change from co2_conc_zero to co2_conc_now
        I.e. starting guess span for bisection
        """

    @abstractmethod
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
