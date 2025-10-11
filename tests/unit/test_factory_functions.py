import pytest

from ciceroscm.carbon_cycle.carbon_cycle_mod import (
    CarbonCycleModel as DefaultCarbonCycleModel,
)
from ciceroscm.carbon_cycle.carbon_cycle_mod_box import (
    CarbonCycleModel as BoxCarbonCycleModel,
)
from ciceroscm.component_factory_functions import (
    create_carbon_cycle_model,
    create_thermal_model,
)
from ciceroscm.thermal_model.two_layer_ocean import TwoLayerOceanModel
from ciceroscm.thermal_model.upwelling_diffusion_model import UpwellingDiffusionModel


def test_create_default_carbon_cycle_model():
    pamset = {"nystart": 1750, "nyend": 2100}
    model = create_carbon_cycle_model("default", pamset)
    assert isinstance(model, DefaultCarbonCycleModel)


def test_create_box_carbon_cycle_model():
    pamset = {"nystart": 1750, "nyend": 2100}
    model = create_carbon_cycle_model("box", pamset)
    assert isinstance(model, BoxCarbonCycleModel)


def test_error_handling_carbon_cycle():
    pamset = {}
    with pytest.raises(ValueError, match="Unknown model type: nonexistant"):
        model = create_carbon_cycle_model("nonexistant", pamset)
        assert model is None


def test_create_default_thermal_model():

    model = create_thermal_model("default")
    assert isinstance(model({}), UpwellingDiffusionModel)


def test_create_twolayer_model():
    model = create_thermal_model("twolayer")
    assert isinstance(model({}), TwoLayerOceanModel)


def test_error_handling_thermal():
    with pytest.raises(ValueError, match="Unknown model type: nonexistant"):
        model = create_thermal_model("nonexistant")
        assert model is None
