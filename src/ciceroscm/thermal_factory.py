"""Thermal model factory function"""

from .two_layer_ocean import TwoLayerOceanModel
from .upwelling_diffusion_model import UpwellingDiffusionModel


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
