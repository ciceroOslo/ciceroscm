"""
Test flat carbon cycle parameter functionality.

Tests the new flat parameter system that allows specifying rb_function and rs_function
parameters as individual float values (e.g., rb_coef0, rb_coef1, rb_tim0, rb_tim1)
instead of requiring dictionary structures.
"""

import pytest

from ciceroscm.carbon_cycle.carbon_cycle_mod import CarbonCycleModel
from ciceroscm.carbon_cycle.rfuns import _process_flat_carbon_parameters


def test_process_flat_rb_parameters():
    """Test processing of flat rb_function parameters."""
    pamset = {
        "rb_coef0": 0.5,
        "rb_coef1": 0.25,
        "rb_coef2": 0.25,
        "rb_tim0": 2.5,
        "rb_tim1": 10.0,
        "rb_tim2": 60.0,
        "other_param": 123.45,
    }

    result = _process_flat_carbon_parameters(pamset)

    # Check rb_function dictionary creation
    assert "rb_function" in result
    assert result["rb_function"]["coeffs"] == [0.5, 0.25, 0.25]
    assert result["rb_function"]["timescales"] == [2.5, 10.0, 60.0]

    # Check flat parameters removed
    for param in ["rb_coef0", "rb_coef1", "rb_coef2", "rb_tim0", "rb_tim1", "rb_tim2"]:
        assert param not in result

    # Check other parameters preserved
    assert result["other_param"] == 123.45


def test_process_flat_rs_parameters():
    """Test processing of flat rs_function parameters."""
    pamset = {
        "rs_coef0": 0.1,
        "rs_coef1": 0.6,
        "rs_coef2": 0.15,
        "rs_coef3": 0.15,
        "rs_tim0": 0.8,
        "rs_tim1": 7.0,
        "rs_tim2": 80.0,
        "other_param": 456.78,
    }

    result = _process_flat_carbon_parameters(pamset)

    # Check rs_function dictionary creation
    assert "rs_function" in result
    assert result["rs_function"]["coeffs"] == [0.1, 0.6, 0.15, 0.15]
    assert result["rs_function"]["timescales"] == [0.8, 7.0, 80.0]

    # Check flat parameters removed
    for param in [
        "rs_coef0",
        "rs_coef1",
        "rs_coef2",
        "rs_coef3",
        "rs_tim0",
        "rs_tim1",
        "rs_tim2",
    ]:
        assert param not in result

    # Check other parameters preserved
    assert result["other_param"] == 456.78


def test_mixed_parameter_styles():
    """Test mixing flat and traditional dictionary parameters."""
    pamset = {
        "rb_coef0": 0.7,
        "rb_coef1": 0.3,
        "rb_tim0": 5.0,
        "rb_tim1": 25.0,
        "rs_function": {"coeffs": [0.2, 0.8], "timescales": [1.0]},
        "other_param": 789.12,
    }

    result = _process_flat_carbon_parameters(pamset)

    # Check rb_function created from flat params
    assert result["rb_function"]["coeffs"] == [0.7, 0.3]
    assert result["rb_function"]["timescales"] == [5.0, 25.0]

    # Check rs_function preserved from dictionary
    assert result["rs_function"]["coeffs"] == [0.2, 0.8]
    assert result["rs_function"]["timescales"] == [1.0]


def test_rb_parameter_mismatch_error():
    """Test error for mismatched rb coefficient and timescale counts."""
    pamset = {
        "rb_coef0": 0.5,
        "rb_coef1": 0.25,
        "rb_tim0": 2.5,
        # Missing rb_tim1
    }

    with pytest.raises(ValueError, match="must match number of rb_tim parameters"):
        _process_flat_carbon_parameters(pamset)


def test_rs_parameter_ratio_error():
    """Test error for incorrect rs coefficient to timescale ratio."""
    pamset = {
        "rs_coef0": 0.1,
        "rs_coef1": 0.6,
        "rs_tim0": 0.8,
        "rs_tim1": 7.0,
        # rs should have one more coef than timescales
    }

    with pytest.raises(
        ValueError, match="must be one more than number of rs_tim parameters"
    ):
        _process_flat_carbon_parameters(pamset)


def test_carbon_cycle_model_integration():
    """Test CarbonCycleModel works with flat parameters."""
    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,  # Short test period
    }

    pamset_carbon_flat = {
        "rb_coef0": 0.5,
        "rb_coef1": 0.25,
        "rb_coef2": 0.25,
        "rb_tim0": 2.5,
        "rb_tim1": 10.0,
        "rb_tim2": 60.0,
        "rs_coef0": 0.1,
        "rs_coef1": 0.6,
        "rs_coef2": 0.15,
        "rs_coef3": 0.15,
        "rs_tim0": 0.8,
        "rs_tim1": 7.0,
        "rs_tim2": 80.0,
        "beta_f": 1.0,
    }

    # Should create model without errors
    model = CarbonCycleModel(pamset_emiconc, pamset_carbon_flat)

    # Check parameters processed correctly
    assert "rb_function" in model.pamset
    assert "rs_function" in model.pamset
    assert model.pamset["rb_function"]["coeffs"] == [0.5, 0.25, 0.25]
    assert model.pamset["rs_function"]["coeffs"] == [0.1, 0.6, 0.15, 0.15]

    # Check flat parameters removed
    assert "rb_coef0" not in model.pamset
    assert "rs_coef0" not in model.pamset


def test_backward_compatibility():
    """Test that traditional dictionary parameters still work."""
    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Traditional dictionary format
    pamset_carbon_dict = {
        "rb_function": {"coeffs": [0.5, 0.25, 0.25], "timescales": [2.5, 10, 60]},
        "rs_function": {"coeffs": [0.1, 0.6, 0.15, 0.15], "timescales": [0.8, 7, 80]},
        "beta_f": 1.0,
    }

    # New flat format
    pamset_carbon_flat = {
        "rb_coef0": 0.5,
        "rb_coef1": 0.25,
        "rb_coef2": 0.25,
        "rb_tim0": 2.5,
        "rb_tim1": 10.0,
        "rb_tim2": 60.0,
        "rs_coef0": 0.1,
        "rs_coef1": 0.6,
        "rs_coef2": 0.15,
        "rs_coef3": 0.15,
        "rs_tim0": 0.8,
        "rs_tim1": 7.0,
        "rs_tim2": 80.0,
        "beta_f": 1.0,
    }

    model_dict = CarbonCycleModel(pamset_emiconc, pamset_carbon_dict)
    model_flat = CarbonCycleModel(pamset_emiconc, pamset_carbon_flat)

    # Both should produce identical results
    assert (
        model_dict.pamset["rb_function"]["coeffs"]
        == model_flat.pamset["rb_function"]["coeffs"]
    )
    assert (
        model_dict.pamset["rb_function"]["timescales"]
        == model_flat.pamset["rb_function"]["timescales"]
    )
    assert (
        model_dict.pamset["rs_function"]["coeffs"]
        == model_flat.pamset["rs_function"]["coeffs"]
    )
    assert (
        model_dict.pamset["rs_function"]["timescales"]
        == model_flat.pamset["rs_function"]["timescales"]
    )


def test_variable_timescale_counts():
    """Test that variable numbers of timescales work correctly."""
    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Test 2-component rb_function
    pamset_2comp = {
        "rb_coef0": 0.7,
        "rb_coef1": 0.3,
        "rb_tim0": 5.0,
        "rb_tim1": 25.0,
        "beta_f": 1.0,
    }

    model_2comp = CarbonCycleModel(pamset_emiconc, pamset_2comp)
    print(model_2comp.pamset)
    assert len(model_2comp.pamset["rb_function"]["coeffs"]) == 2
    assert len(model_2comp.pamset["rb_function"]["timescales"]) == 2

    # Test 4-component rb_function
    pamset_4comp = {
        "rb_coef0": 0.4,
        "rb_coef1": 0.3,
        "rb_coef2": 0.2,
        "rb_coef3": 0.1,
        "rb_tim0": 1.0,
        "rb_tim1": 5.0,
        "rb_tim2": 25.0,
        "rb_tim3": 100.0,
        "beta_f": 1.0,
    }

    model_4comp = CarbonCycleModel(pamset_emiconc, pamset_4comp)
    assert len(model_4comp.pamset["rb_function"]["coeffs"]) == 4
    assert len(model_4comp.pamset["rb_function"]["timescales"]) == 4


def test_rfunction_recomputation_on_reset():
    """Test that r_functions are recomputed when reset_co2_hold is called with new parameters."""
    import numpy as np

    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Initial parameters
    initial_pamset = {
        "rb_coef0": 0.5,
        "rb_coef1": 0.5,
        "rb_tim0": 5.0,
        "rb_tim1": 25.0,
        "beta_f": 1.0,
    }

    model = CarbonCycleModel(pamset_emiconc, initial_pamset)

    # Store initial r_functions
    initial_rb_function = model.r_functions[1, :].copy()

    # Update parameters with different values
    updated_paramset = {
        "rb_coef0": 0.8,  # Changed from 0.5
        "rb_coef1": 0.2,  # Changed from 0.5
        "rb_tim0": 2.0,  # Changed from 5.0
        "rb_tim1": 15.0,  # Changed from 25.0
    }

    # Reset with new parameters
    model.reset_co2_hold(updated_paramset)

    # Get updated r_functions
    updated_rb_function = model.r_functions[1, :].copy()

    # Verify parameter updates were applied
    assert model.pamset["rb_function"]["coeffs"] == [0.8, 0.2]
    assert model.pamset["rb_function"]["timescales"] == [2.0, 15.0]

    # Verify that r_functions actually changed
    assert not np.allclose(
        initial_rb_function, updated_rb_function
    ), "rb_function should have changed with new parameters"


def test_flat_parameters_only_trigger_recomputation():
    """Test that providing only flat parameters triggers r_function recomputation."""
    import numpy as np

    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Create model with no special carbon cycle parameters (uses defaults)
    model = CarbonCycleModel(pamset_emiconc, None)

    # Store initial r_functions (these should be the default functions)
    initial_rb_function = model.r_functions[1, :].copy()
    initial_rs_function = model.r_functions[0, :].copy()

    # Update with ONLY flat parameters (no other carbon cycle parameters)
    flat_params_only = {
        "rb_coef0": 0.8,  # Different from defaults
        "rb_coef1": 0.2,
        "rb_tim0": 3.0,
        "rb_tim1": 15.0,
    }

    # Reset with only flat parameters
    model.reset_co2_hold(flat_params_only)

    # Get updated r_functions
    updated_rb_function = model.r_functions[1, :].copy()
    updated_rs_function = model.r_functions[0, :].copy()

    # Check that rb_function changed (because we provided new parameters)
    rb_functions_changed = not np.allclose(initial_rb_function, updated_rb_function)

    # Check that rs_function did NOT change (because we didn't provide rs parameters)
    rs_functions_unchanged = np.allclose(initial_rs_function, updated_rs_function)

    # Verify parameter updates were applied
    assert "rb_function" in model.pamset
    print(model.pamset["rb_function"])
    assert model.pamset["rb_function"]["coeffs"] == [0.8, 0.2]
    assert model.pamset["rb_function"]["timescales"] == [3.0, 15.0]

    # Verify that rb_function changed but rs_function didn't
    assert rb_functions_changed, "rb_function should have changed with new parameters"
    assert rs_functions_unchanged, "rs_function should not have changed"


def test_non_function_parameters_dont_trigger_recomputation():
    """Test that non-function parameters don't trigger unnecessary r_function recomputation."""
    import numpy as np

    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    model = CarbonCycleModel(pamset_emiconc, None)

    initial_rb_function = model.r_functions[1, :].copy()
    initial_rs_function = model.r_functions[0, :].copy()

    # Update with only non-function parameters
    non_function_params = {
        "beta_f": 0.5,  # Changed from default
        "mixed_carbon": 70,  # Changed from default
    }

    model.reset_co2_hold(non_function_params)

    updated_rb_function = model.r_functions[1, :].copy()
    updated_rs_function = model.r_functions[0, :].copy()

    # Check that neither function changed
    rb_functions_unchanged = np.allclose(initial_rb_function, updated_rb_function)
    rs_functions_unchanged = np.allclose(initial_rs_function, updated_rs_function)

    # Verify parameter updates were applied
    assert model.pamset["beta_f"] == 0.5
    assert model.pamset["mixed_carbon"] == 70

    # Verify that neither function changed
    assert rb_functions_unchanged, "rb_function should not have changed"
    assert rs_functions_unchanged, "rs_function should not have changed"
