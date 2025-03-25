"""Carbon cycle factory function"""

from .carbon_cycle_mod import CarbonCycleModel as DefaultCarbonCycleModel
from .carbon_cycle_mod_box import CarbonCycleModel as BoxCarbonCycleModel


def create_carbon_cycle_model(model_type, pamset):
    """
    Create a CarbonCycleModel instance (factory function).

    Parameters
    ----------
    model_type : str
        The type of carbon cycle model to create ("default" or "box").
    pamset : dict
        Parameter set for the model.

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
        return DefaultCarbonCycleModel(pamset)
    if model_type == "box":
        return BoxCarbonCycleModel(pamset)
    raise ValueError(f"Unknown model type: {model_type}")
