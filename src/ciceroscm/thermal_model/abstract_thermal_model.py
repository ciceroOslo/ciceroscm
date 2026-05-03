"""
Abstract carbon cycle model
"""

from abc import ABC, abstractmethod

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

    @classmethod
    def get_thermal_model_required_pamset(cls):
        """
        Get the class variable carbon_cycle_model_required_pamset
        which should be a class variable. This setup is to
        substitute for the lack of abstract class variables in python
        """
        return cls.thermal_model_required_pamset

    @property
    @abstractmethod
    def output_dict_default(self):
        """
        Overwrite this to set a default output dictionary
        structure for the thermal model
        """

    @classmethod
    def get_output_dict_thermal(cls):
        """
        Get the class variable carbon_cycle_model_required_pamset
        which should be a class variable. This setup is to
        substitute for the lack of abstract class variables in python
        """
        return cls.output_dict_default

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
            self.get_thermal_model_required_pamset(),
            pamset,
            cut_warnings=True,
        )

    # TODO: Should this have an option to pass the state?
    # Alternatively, should there be a separate method to
    # update the state?
    # Also now we send four specific forcings, maybe it should be a dict?
    @abstractmethod
    def energy_budget(self, forc_nh, forc_sh, fn_volc, fs_volc):
        """
        Calculate the energy balance

        Parameters
        ----------
        forc_nh : float
            Northern hemispheric forcing (W/m^2)
        forc_sh : float
            Southern hemispheric forcing (W/m^2)
        fn_volc : float
            Northern hemispheric volcanic forcing (W/m^2)
        fs_volc : float
            Northern hemispheric volcanic forcing (W/m^2)

        Returns
        -------
        float
            The energy balance
        """

    # ------------------------------------------------------------------
    # Optional pattern-mediated feedback interface.
    #
    # Thermal models that support a forcing-composition-dependent feedback
    # parameter override these two methods. The driver only invokes them
    # when delta_lambda_aero != 0; otherwise existing behaviour is
    # unchanged. The interface is intentionally in Gregory units
    # (W m^-2 K^-1) so the driver does not need to know each model's
    # internal sign/inverse convention for storing lambda.
    # ------------------------------------------------------------------
    def get_feedback_gregory(self):
        """
        Return the current climate feedback parameter in Gregory units.

        Returned value should be a float.
        Feedback parameter (W m^-2 K^-1), positive for a stabilising
        feedback. Override in concrete thermal models that support
        pattern-mediated feedback modulation.
        """
        return None

    def set_feedback_gregory(self, w_aero):
        """
        Update the climate feedback parameter to ``lambda_eff``.

        Used by the driver to apply pattern-mediated feedback updates
        each year (Tier 3 of the pattern-effect formulation). Implementations
        must update any cached derived quantities that depend on lambda
        so that the next ``energy_budget`` call uses the new value.

        Parameters
        ----------
        w_aero : float
            New aerosol based unitless feedback parameter
            To be multiplied with delta_lambda_aero which is in Gregory units (W m^-2 K^-1)
            And added to the ECS parameter.

        """
