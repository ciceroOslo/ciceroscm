"""Factory functions, currently for carbon cycle and thermal model"""

from .carbon_cycle.carbon_cycle_mod import CarbonCycleModel as DefaultCarbonCycleModel
from .carbon_cycle.carbon_cycle_mod_box import CarbonCycleModel as BoxCarbonCycleModel
from .thermal_model.two_layer_ocean import TwoLayerOceanModel
from .thermal_model.upwelling_diffusion_model import UpwellingDiffusionModel


def create_carbon_cycle_model(model_type, pamset, pamset_carbon=None):
    """
    Create a CarbonCycleModel instance (factory function).

    Parameters
    ----------
    model_type : str
        The type of carbon cycle model to create ("default" or "box").
    pamset : dict
        Parameter set for the model, shared with concentrations emissions
        handler
    pamset_carbon : dict
        Carbon cycle model parameters

    Returns
    -------
    CarbonCycleModel
        An instance of the selected carbon cycle model.

    Raises
    ------
    ValueError
        If an undefined carbon cycle model_type is called for
    """
    if model_type == "default":
        return DefaultCarbonCycleModel(pamset, pamset_carbon)
    if model_type == "box":
        return BoxCarbonCycleModel(pamset, pamset_carbon)
    raise ValueError(f"Unknown model type: {model_type}")


def create_thermal_model(model_type):
    """
    Create a ThermalModel instance (factory function).

    Parameters
    ----------
    model_type : str
        The type of carbon cycle model to create ("default" or "box").

    Returns
    -------
    CarbonCycleModel
        An instance of the selected carbon cycle model.

    Raises
    ------
    ValueError
        If an undefined thermal model_type is called for
    """
    if model_type == "default":
        return UpwellingDiffusionModel
    if model_type == "twolayer":
        return TwoLayerOceanModel
    raise ValueError(f"Unknown model type: {model_type}")
