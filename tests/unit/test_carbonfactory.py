import pytest

from ciceroscm.carbon_cycle_factory import create_carbon_cycle_model
from ciceroscm.carbon_cycle_mod import CarbonCycleModel as DefaultCarbonCycleModel
from ciceroscm.carbon_cycle_mod_box import CarbonCycleModel as BoxCarbonCycleModel


def test_create_default_carbon_cycle_model():
    pamset = {"nystart": 1750, "nyend": 2100}
    model = create_carbon_cycle_model("default", pamset)
    assert isinstance(model, DefaultCarbonCycleModel)


def test_create_box_carbon_cycle_model():
    pamset = {"nystart": 1750, "nyend": 2100}
    model = create_carbon_cycle_model("box", pamset)
    assert isinstance(model, BoxCarbonCycleModel)


def test_error_handling():
    pamset = {}
    with pytest.raises(ValueError, match="Unknown model type: nonexistant"):
        model = create_carbon_cycle_model("nonexistant", pamset)
        assert model is None
