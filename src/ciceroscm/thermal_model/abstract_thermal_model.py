

"""
Abstract carbon cycle model
"""

from abc import ABC, abstractmethod

import numpy as np

from .._utils import cut_and_check_pamset



class AbstractThermalModel(ABC):
    """
    Abstract class to define carbon cycle methodology
    """

    @property
    @abstractmethod
    def thermal_model_required_pamset(self):
        """
        Overwrite this to set a required thermal
        model parameterset dictionary. This should be
        a class variable
        """
        pass

    @classmethod
    def get_thermal_model_required_pamset(cls):
        """
        Get the class variable carbon_cycle_model_required_pamset
        which should be a class variable. This setup is to
        substitute for the lack of abstract class variables in python
        """
        return cls.thermal_model_required_pamset

    def __init__(self, pamset=None):
        """
        Initialize the thermal model

        Parameters
        ----------
        pamset : dict, optional
            A dictionary of parameters for the thermal model.
            The default is None.
        """
        if pamset is None:
            pamset = {}

        self.pamset = cut_and_check_pamset(
            pamset,
            self.get_thermal_model_required_pamset(),
            cut_warnings=True,
        )

    # TODO: Should this have an option to pass the state?
    # Alternatively, should there be a separate method to
    # update the state?
    @abstractmethod
    def energy_balance(self, forcings):
        """
        Calculate the energy balance

        Parameters
        ----------
        forcings : dict
            A dictionary of forcings

        Returns
        -------
        float
            The energy balance
        """
        pass