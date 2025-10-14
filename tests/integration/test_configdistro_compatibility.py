"""
Test _ConfigDistro compatibility with flat carbon cycle parameters.

Tests if flat rb_/rs_ parameters are properly handled by _ConfigDistro
"""

import pytest

from ciceroscm.carbon_cycle.carbon_cycle_mod import CarbonCycleModel
from ciceroscm.parallel._configdistro import _ConfigDistro, prior_flat


def test_configdistro_flat_parameters_phase1_fix():
    """
    Test that configdistro flat carbon parameters to pamset_carbon.

    Requirements:
    1. Flat parameters are included in default prior_flat
    2. rb_/rs_ parameters are explicitly routed to pamset_carbon
    3. _process_flat_carbon_parameters correctly processes them
    4. CarbonCycleModel is created with proper rb_function/rs_function
    """
    # Test with default _ConfigDistro (now includes flat parameters)
    config_distro = _ConfigDistro()  # Uses default prior_flat with flat parameters
    config_list = config_distro.make_config_lists(1)
    config = config_list[0]

    # Verify flat parameters now correctly go to pamset_carbon
    flat_params_emiconc = {
        k: v
        for k, v in config["pamset_emiconc"].items()
        if k.startswith(("rb_", "rs_"))
    }
    flat_params_carbon = {
        k: v for k, v in config["pamset_carbon"].items() if k.startswith(("rb_", "rs_"))
    }

    assert (
        len(flat_params_carbon) > 0
    ), "Flat parameters should be in pamset_carbon (fixed behavior)"
    assert (
        len(flat_params_emiconc) == 0
    ), "Flat parameters should NOT be in pamset_emiconc (fixed behavior)"

    # Try to create CarbonCycleModel as it would normally be used
    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
        **config["pamset_emiconc"],
    }

    model = CarbonCycleModel(pamset_emiconc, config["pamset_carbon"])

    # Fixed behavior: rb_function and rs_function ARE created
    assert (
        "rb_function" in model.pamset
    ), "rb_function should be created when flat params are in pamset_carbon"
    assert (
        "rs_function" in model.pamset
    ), "rs_function should be created when flat params are in pamset_carbon"

    # Verify functions are properly constructed
    rb_func = model.pamset["rb_function"]
    rs_func = model.pamset["rs_function"]

    assert (
        "coeffs" in rb_func and "timescales" in rb_func
    ), "rb_function should have coeffs and timescales"
    assert (
        "coeffs" in rs_func and "timescales" in rs_func
    ), "rs_function should have coeffs and timescales"
    assert len(rb_func["coeffs"]) == len(
        rb_func["timescales"]
    ), "rb_function coeffs and timescales should match"
    assert (
        len(rs_func["coeffs"]) == len(rs_func["timescales"]) + 1
    ), "rs_function should have one more coeff than timescales"


def test_configdistro_flat_parameters_with_custom_params():
    """
    Test that custom flat parameters are also correctly routed.

    This test verifies that when users provide custom rb_/rs_ parameters via
    distro_dict, they are correctly routed to pamset_carbon and processed.
    """
    # Create custom distribution with flat carbon parameters
    custom_distro = prior_flat.copy()
    custom_distro.update(
        {
            "rb_coef0": [0.3, 0.7],
            "rb_tim0": [1.0, 5.0],
            "rs_coef0": [0.3, 0.6],
        }
    )

    config_distro = _ConfigDistro(distro_dict=custom_distro)
    config_list = config_distro.make_config_lists(1)
    config = config_list[0]

    #  custom flat parameters go to pamset_carbon
    flat_params_emiconc = {
        k: v
        for k, v in config["pamset_emiconc"].items()
        if k.startswith(("rb_", "rs_"))
    }
    flat_params_carbon = {
        k: v for k, v in config["pamset_carbon"].items() if k.startswith(("rb_", "rs_"))
    }
    print(flat_params_carbon)
    assert len(flat_params_emiconc) == 0, "No flat params should be in pamset_emiconc"
    assert len(flat_params_carbon) == 13, "Flat params should be in pamset_carbon"

    # Create model directly
    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
        **config["pamset_emiconc"],
    }

    model = CarbonCycleModel(pamset_emiconc, config["pamset_carbon"])

    # rb_function should be created correctly
    assert "rb_function" in model.pamset, "rb_function should be created"

    # Since we specified rb_coef0, rb_coef1, rb_tim0, rb_tim1, but the default
    # prior_flat also includes rb_coef2, rb_tim2, we should get 3 coefficients/timescales
    rb_func = model.pamset["rb_function"]
    assert len(rb_func["coeffs"]) == 3, "Should have 3 rb coefficients (rb_coef0-2)"
    assert len(rb_func["timescales"]) == 3, "Should have 3 rb timescales (rb_tim0-2)"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
