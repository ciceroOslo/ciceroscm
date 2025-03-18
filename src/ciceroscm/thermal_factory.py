from .upwelling_diffusion_model import UpwellingDiffusionModel
from .two_layer_ocean import TwoLayerOceanModel


def create_thermal_model(model_type):
    """
    Factory function to create a ThermalModel instance.

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
        return UpwellingDiffusionModel
    elif model_type == "twolayer":
        return TwoLayerOceanModel
    else:
        raise ValueError(f"Unknown model type: {model_type}")
