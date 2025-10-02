#!/usr/bin/env python3
"""
Test to verify that flat parameters are properly recognized in precalc_r_functions.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

import numpy as np
from ciceroscm.carbon_cycle.carbon_cycle_mod import CarbonCycleModel


def test_initialization_with_flat_params():
    """Test that initialization with flat parameters correctly uses custom functions."""
    print("Testing initialization with flat parameters...")

    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Create model with only flat rs parameters
    pamset_with_flat_rs = {
        "rs_coef0": 0.2,
        "rs_coef1": 0.5,
        "rs_coef2": 0.3,
        "rs_tim0": 1.5,
        "rs_tim1": 8.0,
        "beta_f": 1.0,
    }

    model = CarbonCycleModel(pamset_emiconc, pamset_with_flat_rs)

    # Check that rs_function was properly added to pamset
    assert (
        "rs_function" in model.pamset
    ), "rs_function should be in pamset after flat parameter processing"
    assert model.pamset["rs_function"]["coeffs"] == [0.2, 0.5, 0.3]
    assert model.pamset["rs_function"]["timescales"] == [1.5, 8.0]

    # Create a model with default parameters to compare
    model_default = CarbonCycleModel(pamset_emiconc, None)

    # The r_functions should be different
    rs_functions_different = not np.allclose(
        model.r_functions[0, :5], model_default.r_functions[0, :5]
    )

    print(f"Custom rs_function[0:5]: {model.r_functions[0, :5]}")
    print(f"Default rs_function[0:5]: {model_default.r_functions[0, :5]}")
    print(f"Functions are different: {rs_functions_different}")

    assert rs_functions_different, "Custom rs_function should be different from default"

    print("✓ Initialization with flat parameters works correctly")


def test_initialization_with_flat_rb_params():
    """Test that initialization with flat rb parameters correctly uses custom functions."""
    print("\nTesting initialization with flat rb parameters...")

    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Create model with only flat rb parameters
    pamset_with_flat_rb = {
        "rb_coef0": 0.8,
        "rb_coef1": 0.2,
        "rb_tim0": 3.0,
        "rb_tim1": 15.0,
        "beta_f": 1.0,
    }

    model = CarbonCycleModel(pamset_emiconc, pamset_with_flat_rb)

    # Check that rb_function was properly added to pamset
    assert (
        "rb_function" in model.pamset
    ), "rb_function should be in pamset after flat parameter processing"
    assert model.pamset["rb_function"]["coeffs"] == [0.8, 0.2]
    assert model.pamset["rb_function"]["timescales"] == [3.0, 15.0]

    # Create a model with default parameters to compare
    model_default = CarbonCycleModel(pamset_emiconc, None)

    # The r_functions should be different
    rb_functions_different = not np.allclose(
        model.r_functions[1, :5], model_default.r_functions[1, :5]
    )

    print(f"Custom rb_function[0:5]: {model.r_functions[1, :5]}")
    print(f"Default rb_function[0:5]: {model_default.r_functions[1, :5]}")
    print(f"Functions are different: {rb_functions_different}")

    assert rb_functions_different, "Custom rb_function should be different from default"

    print("✓ Initialization with flat rb parameters works correctly")


def test_runtime_update_recognition():
    """Test that runtime updates properly recognize when to use custom vs default functions."""
    print("\nTesting runtime update recognition...")

    pamset_emiconc = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }

    # Start with a model that has no custom functions (uses defaults)
    model = CarbonCycleModel(pamset_emiconc, None)

    # Check that it's using default functions (no rs_function or rb_function in pamset)
    assert "rs_function" not in model.pamset, "Should start with no custom rs_function"
    assert "rb_function" not in model.pamset, "Should start with no custom rb_function"

    # Store initial functions
    initial_rs = model.r_functions[0, :5].copy()
    initial_rb = model.r_functions[1, :5].copy()

    print(f"Initial rs_function[0:5]: {initial_rs}")
    print(f"Initial rb_function[0:5]: {initial_rb}")

    # Update with flat parameters - this should add rs_function to pamset and trigger recomputation
    model.reset_co2_hold(
        {
            "rs_coef0": 0.3,
            "rs_coef1": 0.7,
            "rs_tim0": 2.0,
        }
    )

    # Check that rs_function was added to pamset
    assert "rs_function" in model.pamset, "rs_function should be added to pamset"
    assert model.pamset["rs_function"]["coeffs"] == [0.3, 0.7]
    assert model.pamset["rs_function"]["timescales"] == [2.0]

    # Check that rb_function is still not in pamset (we didn't update it)
    assert (
        "rb_function" not in model.pamset
    ), "rb_function should still not be in pamset"

    # Check the functions changed appropriately
    updated_rs = model.r_functions[0, :5].copy()
    updated_rb = model.r_functions[1, :5].copy()

    print(f"Updated rs_function[0:5]: {updated_rs}")
    print(f"Updated rb_function[0:5]: {updated_rb}")

    rs_changed = not np.allclose(initial_rs, updated_rs)
    rb_unchanged = np.allclose(initial_rb, updated_rb)

    print(f"rs_function changed: {rs_changed}")
    print(f"rb_function unchanged: {rb_unchanged}")

    assert rs_changed, "rs_function should have changed"
    assert rb_unchanged, "rb_function should not have changed"

    print("✓ Runtime update recognition works correctly")


def main():
    """Run all tests."""
    print("Testing flat parameter recognition in precalc_r_functions...\n")

    test_initialization_with_flat_params()
    test_initialization_with_flat_rb_params()
    test_runtime_update_recognition()

    print("\n🎉 All flat parameter recognition tests passed!")


if __name__ == "__main__":
    main()
