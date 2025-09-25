from abc import ABC, abstractmethod

from .._utils import update_pam_if_numeric
from common_carbon_cycle_functions import carbon_cycle_init_pamsets, calculate_airborne_fraction


class AbstractCarbonCycleModel(ABC):
    """
    Abstract class to define carbon cycle methodology
    """
    def carbon_cycle_model_required_pamset(self):
        pass

    def get_carbon_cycle_required_pamset(cls):
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
            pamset_emiconc,
            pamset_carbon,
            self.get_carbon_cycle_required_pamset()
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
                can_change=self.get_carbon_cycle_required_pamset.keys(),
            )

    # TODO: Generalise so more than just temperature can be sent from thermal
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
        pass

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
        pass
