"""
Test that flat carbon parameters are actually perturbed in distribution runs
and that they affect model output.

This test demonstrates that:
1. Flat parameters vary across different samples in a distribution run
2. Different flat parameter values produce different model outputs
3. The carbon cycle responds to flat parameter perturbations
"""

import pytest
import numpy as np

from ciceroscm.parallel._configdistro import _ConfigDistro, prior_flat
from ciceroscm.carbon_cycle.carbon_cycle_mod import CarbonCycleModel, PREINDUSTRIAL_CO2_CONC


def test_flat_parameters_are_perturbed_in_distro():
    """Test that flat parameters actually vary across multiple samples."""
    # Create multiple samples with flat parameters
    config_distro = _ConfigDistro()  # Uses default prior_flat with flat parameters
    n_samples = 10
    config_list = config_distro.make_config_lists(n_samples)
    
    # Extract flat parameter values from all samples
    all_rb_coef0 = []
    all_rs_coef0 = []
    
    for config in config_list:
        pamset_carbon = config['pamset_carbon']
        if 'rb_coef0' in pamset_carbon:
            all_rb_coef0.append(pamset_carbon['rb_coef0'])
        if 'rs_coef0' in pamset_carbon:
            all_rs_coef0.append(pamset_carbon['rs_coef0'])
    
    # Verify we have flat parameters in all samples
    assert len(all_rb_coef0) == n_samples, "rb_coef0 should be present in all samples"
    assert len(all_rs_coef0) == n_samples, "rs_coef0 should be present in all samples"
    
    # Verify flat parameters actually vary across samples
    rb_coef0_std = np.std(all_rb_coef0)
    rs_coef0_std = np.std(all_rs_coef0)
    
    assert rb_coef0_std > 0, f"rb_coef0 should vary across samples, got std={rb_coef0_std}"
    assert rs_coef0_std > 0, f"rs_coef0 should vary across samples, got std={rs_coef0_std}"
    
    # Verify values are within expected ranges from prior_flat
    rb_coef0_min, rb_coef0_max = prior_flat['rb_coef0']
    rs_coef0_min, rs_coef0_max = prior_flat['rs_coef0']
    
    assert all(rb_coef0_min <= val <= rb_coef0_max for val in all_rb_coef0), \
        f"rb_coef0 values should be within [{rb_coef0_min}, {rb_coef0_max}]"
    assert all(rs_coef0_min <= val <= rs_coef0_max for val in all_rs_coef0), \
        f"rs_coef0 values should be within [{rs_coef0_min}, {rs_coef0_max}]"
    
    print(f"✅ rb_coef0 varies: {min(all_rb_coef0):.3f} to {max(all_rb_coef0):.3f} (std={rb_coef0_std:.3f})")
    print(f"✅ rs_coef0 varies: {min(all_rs_coef0):.3f} to {max(all_rs_coef0):.3f} (std={rs_coef0_std:.3f})")


def test_flat_parameters_affect_carbon_functions():
    """Test that different flat parameters produce different carbon functions."""
    # Create two different parameter sets with contrasting flat parameter values
    config_distro = _ConfigDistro()
    config_list = config_distro.make_config_lists(2)
    
    # Get two different configurations
    config1 = config_list[0]
    config2 = config_list[1]
    
    # Create CarbonCycleModels with the different configurations
    pamset_emiconc_base = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1755,
    }
    
    pamset_emiconc1 = {**pamset_emiconc_base, **config1['pamset_emiconc']}
    pamset_emiconc2 = {**pamset_emiconc_base, **config2['pamset_emiconc']}
    
    model1 = CarbonCycleModel(pamset_emiconc1, config1['pamset_carbon'])
    model2 = CarbonCycleModel(pamset_emiconc2, config2['pamset_carbon'])
    
    # Both should have rb_function and rs_function
    assert 'rb_function' in model1.pamset, "Model 1 should have rb_function"
    assert 'rb_function' in model2.pamset, "Model 2 should have rb_function"
    assert 'rs_function' in model1.pamset, "Model 1 should have rs_function"
    assert 'rs_function' in model2.pamset, "Model 2 should have rs_function"
    
    # Extract the functions
    rb_func1 = model1.pamset['rb_function']
    rb_func2 = model2.pamset['rb_function']
    rs_func1 = model1.pamset['rs_function']
    rs_func2 = model2.pamset['rs_function']
    
    # Functions should be different (since flat parameters are different)
    rb_coeffs_different = not np.allclose(rb_func1['coeffs'], rb_func2['coeffs'])
    rb_times_different = not np.allclose(rb_func1['timescales'], rb_func2['timescales'])
    rs_coeffs_different = not np.allclose(rs_func1['coeffs'], rs_func2['coeffs'])
    rs_times_different = not np.allclose(rs_func1['timescales'], rs_func2['timescales'])
    
    assert rb_coeffs_different or rb_times_different, \
        "rb_function should differ between models with different flat parameters"
    assert rs_coeffs_different or rs_times_different, \
        "rs_function should differ between models with different flat parameters"
    
    print(f"✅ rb_function differs between models:")
    print(f"   Model 1 coeffs: {[f'{c:.3f}' for c in rb_func1['coeffs']]}")
    print(f"   Model 2 coeffs: {[f'{c:.3f}' for c in rb_func2['coeffs']]}")
    print(f"✅ rs_function differs between models:")
    print(f"   Model 1 coeffs: {[f'{c:.3f}' for c in rs_func1['coeffs']]}")
    print(f"   Model 2 coeffs: {[f'{c:.3f}' for c in rs_func2['coeffs']]}")


def test_flat_parameters_affect_model_output():
    """Test that flat parameter perturbations actually affect carbon cycle model output."""
    # Create two configurations with manually controlled flat parameters
    # to ensure they're different enough to see model output differences
    
    # Configuration 1: Fast carbon cycle (smaller timescales, different coeffs)
    custom_distro1 = prior_flat.copy()
    custom_distro1.update({
        "rb_coef0": [0.59, 0.61],  # Narrow range around 0.6 - high first coefficient
        "rb_coef1": [0.19, 0.21],  # Narrow range around 0.2
        "rb_coef2": [0.19, 0.21],  # Narrow range around 0.2
        "rb_tim0": [1.4, 1.6],     # Narrow range around 1.5 - fast timescale
        "rb_tim1": [7.9, 8.1],     # Narrow range around 8.0
        "rb_tim2": [44.0, 46.0],   # Narrow range around 45.0
        "rs_coef0": [0.09, 0.11],  # Narrow range around 0.1
        "rs_coef1": [0.59, 0.61],  # Narrow range around 0.6
        "rs_coef2": [0.14, 0.16],  # Narrow range around 0.15
        "rs_coef3": [0.14, 0.16],  # Narrow range around 0.15
        "rs_tim0": [0.7, 0.9],     # Narrow range around 0.8
        "rs_tim1": [5.9, 6.1],     # Narrow range around 6.0
        "rs_tim2": [64.0, 66.0],   # Narrow range around 65.0
    })
    
    # Configuration 2: Slow carbon cycle (larger timescales, different coeffs)
    custom_distro2 = prior_flat.copy()
    custom_distro2.update({
        "rb_coef0": [0.29, 0.31],  # Narrow range around 0.3 - low first coefficient  
        "rb_coef1": [0.34, 0.36], # Narrow range around 0.35
        "rb_coef2": [0.34, 0.36], # Narrow range around 0.35
        "rb_tim0": [3.9, 4.1],    # Narrow range around 4.0 - slow timescale
        "rb_tim1": [11.9, 12.1],  # Narrow range around 12.0
        "rb_tim2": [69.0, 71.0],  # Narrow range around 70.0
        "rs_coef0": [0.14, 0.16], # Narrow range around 0.15
        "rs_coef1": [0.49, 0.51], # Narrow range around 0.5
        "rs_coef2": [0.17, 0.18], # Narrow range around 0.175
        "rs_coef3": [0.17, 0.18], # Narrow range around 0.175
        "rs_tim0": [1.0, 1.2],    # Narrow range around 1.1
        "rs_tim1": [8.9, 9.1],    # Narrow range around 9.0
        "rs_tim2": [84.0, 86.0],  # Narrow range around 85.0
    })
    
    # Create models with these contrasting parameter sets
    config_distro1 = _ConfigDistro(distro_dict=custom_distro1)
    config_distro2 = _ConfigDistro(distro_dict=custom_distro2)
    
    config1 = config_distro1.make_config_lists(1)[0]
    config2 = config_distro2.make_config_lists(1)[0]
    
    # Base parameter set for both models
    pamset_emiconc_base = {
        "idtm": 24,
        "nystart": 1750,
        "nyend": 1760,  # Run for 10 years to see differences
    }
    
    pamset_emiconc1 = {**pamset_emiconc_base, **config1['pamset_emiconc']}
    pamset_emiconc2 = {**pamset_emiconc_base, **config2['pamset_emiconc']}
    
    model1 = CarbonCycleModel(pamset_emiconc1, config1['pamset_carbon'])
    model2 = CarbonCycleModel(pamset_emiconc2, config2['pamset_carbon'])
    
    # Run both models with the same CO2 emission pulse for a single year
    year = 1750
    emission = 100.0  # 100 GtC pulse
    
    # Get atmospheric CO2 response for both models
    co2_atm1 = model1.co2em2conc(year, emission)
    co2_atm2 = model2.co2em2conc(year, emission)
    
    # The models should produce different atmospheric CO2 concentrations
    # due to different carbon cycle parameters
    abs_diff = abs(co2_atm1 - co2_atm2)
    
    assert abs_diff > 0.01, \
        f"Models with different flat parameters should produce different outputs, diff={abs_diff:.6f} ppm"
    
    # Check that both models produced reasonable responses
    assert co2_atm1 > PREINDUSTRIAL_CO2_CONC, "Model 1 should show CO2 increase from emission pulse"
    assert co2_atm2 > PREINDUSTRIAL_CO2_CONC, "Model 2 should show CO2 increase from emission pulse"
    
    print(f"✅ Model outputs differ significantly:")
    print(f"   Model 1 CO2 concentration: {co2_atm1:.3f} ppm")
    print(f"   Model 2 CO2 concentration: {co2_atm2:.3f} ppm")
    print(f"   Absolute difference: {abs_diff:.3f} ppm")
    
    # Print the actual flat parameter values that caused the difference
    print(f"✅ Flat parameter values that caused the difference:")
    for param in ['rb_coef0', 'rb_tim0', 'rs_coef0', 'rs_tim0']:
        val1 = config1['pamset_carbon'].get(param, 'N/A')
        val2 = config2['pamset_carbon'].get(param, 'N/A')
        print(f"   {param}: Model1={val1:.3f}, Model2={val2:.3f}")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])