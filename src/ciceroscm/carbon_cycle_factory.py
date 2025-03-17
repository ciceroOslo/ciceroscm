from .carbon_cycle_mod import CarbonCycleModel as DefaultCarbonCycleModel
from .carbon_cycle_mod_box import CarbonCycleModel as BoxCarbonCycleModel

def create_carbon_cycle_model(model_type, pamset):
    """
    Factory function to create a CarbonCycleModel instance.

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
    """
    if model_type == "default":
        return DefaultCarbonCycleModel(pamset)
    elif model_type == "box":
        return BoxCarbonCycleModel(pamset)
    else:
        raise ValueError(f"Unknown model type: {model_type}")