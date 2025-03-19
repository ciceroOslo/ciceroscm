import pytest

from ciceroscm.thermal_factory import create_thermal_model
from ciceroscm.two_layer_ocean import TwoLayerOceanModel
from ciceroscm.upwelling_diffusion_model import UpwellingDiffusionModel


def test_create_default_thermal_model():

    model = create_thermal_model("default")
    assert isinstance(model({}), UpwellingDiffusionModel)


def test_create_twolayer_model():
    model = create_thermal_model("twolayer")
    assert isinstance(model({}), TwoLayerOceanModel)


def test_error_handling():
    with pytest.raises(ValueError, match="Unknown model type: nonexistant"):
        model = create_thermal_model("nonexistant")
        assert model is None
