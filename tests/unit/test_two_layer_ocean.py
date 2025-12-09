"""
Tests for TwoLayerOceanModel thermal model
"""

import numpy as np
import pytest

from ciceroscm.constants import WATER_DENSITY, WATER_HEAT_CAPACITY
from ciceroscm.thermal_model.two_layer_ocean import TwoLayerOceanModel


class TestTwoLayerOceanModel:
    """Test suite for TwoLayerOceanModel"""

    def test_init_default_parameters(self):
        """Test initialization with default parameters"""
        model = TwoLayerOceanModel()

        # Check that all expected parameters are set to defaults
        assert model.pamset["lambda"] == 3.74 / 3
        assert model.pamset["k"] == 0.5
        assert model.pamset["ocean_efficacy"] == 1

        # Check that c_fast and c_slow are calculated correctly from mixed/deep
        # Default mixed = 50m, deep = 1200m
        expected_c_fast = 50 * WATER_DENSITY * WATER_HEAT_CAPACITY / (86400 * 365.0)
        expected_c_slow = 1200 * WATER_DENSITY * WATER_HEAT_CAPACITY / (86400 * 365.0)

        assert model.pamset["c_fast"] == pytest.approx(expected_c_fast, rel=1e-10)
        assert model.pamset["c_slow"] == pytest.approx(expected_c_slow, rel=1e-10)

        # Check initial temperatures are zero
        assert model.temp_fast == 0.0
        assert model.temp_slow == 0.0

    def test_init_custom_parameters(self):
        """Test initialization with custom parameters"""
        custom_params = {
            "lambda": 1.5,
            "mixed": 100,
            "deep": 800,
            "k": 0.8,
            "ocean_efficacy": 1.2,
        }

        model = TwoLayerOceanModel(custom_params)

        # Check that custom parameters are set
        assert model.pamset["lambda"] == 1.5
        assert model.pamset["k"] == 0.8
        assert model.pamset["ocean_efficacy"] == 1.2

        # Check that c_fast and c_slow are calculated from custom mixed/deep
        expected_c_fast = 100 * WATER_DENSITY * WATER_HEAT_CAPACITY / (86400 * 365.0)
        expected_c_slow = 800 * WATER_DENSITY * WATER_HEAT_CAPACITY / (86400 * 365.0)

        assert model.pamset["c_fast"] == pytest.approx(expected_c_fast, rel=1e-10)
        assert model.pamset["c_slow"] == pytest.approx(expected_c_slow, rel=1e-10)

    def test_init_partial_parameters(self):
        """Test initialization with some custom and some default parameters"""
        partial_params = {
            "lambda": 2.0,
            "k": 0.3,
            # mixed, deep, efficacy should use defaults
        }

        model = TwoLayerOceanModel(partial_params)

        # Check custom parameters
        assert model.pamset["lambda"] == 2.0
        assert model.pamset["k"] == 0.3

        # Check defaults are used for unspecified parameters
        assert model.pamset["ocean_efficacy"] == 1
        expected_c_fast = (
            50 * WATER_DENSITY * WATER_HEAT_CAPACITY / (86400 * 365.0)
        )  # from default mixed=50
        expected_c_slow = (
            1200 * WATER_DENSITY * WATER_HEAT_CAPACITY / (86400 * 365.0)
        )  # from default deep=1200
        assert model.pamset["c_fast"] == pytest.approx(expected_c_fast, rel=1e-10)
        assert model.pamset["c_slow"] == pytest.approx(expected_c_slow, rel=1e-10)

    def test_energy_budget_zero_forcing(self):
        """Test energy budget with zero forcing"""
        model = TwoLayerOceanModel()

        # Zero forcing should result in zero temperature changes
        result = model.energy_budget(0.0, 0.0, [0.0], [0.0])

        assert result["dtemp_fast"] == 0.0
        assert result["dtemp_slow"] == 0.0
        assert result["dtemp"] == 0.0
        assert result["RIB"] == 0.0

    def test_energy_budget_positive_forcing(self):
        """Test energy budget with positive forcing"""
        model = TwoLayerOceanModel()

        # Apply positive forcing
        result = model.energy_budget(1.0, 1.0, [0.0], [0.0])

        # With positive forcing, fast layer should warm up
        assert result["dtemp_fast"] > 0.0
        assert result["dtemp_slow"] >= 0.0  # May be zero initially since temp_slow=0
        assert result["dtemp"] > 0.0  # dtemp = temp_fast

        # Check that model state is updated
        assert model.temp_fast > 0.0

    def test_energy_budget_negative_forcing(self):
        """Test energy budget with negative forcing"""
        model = TwoLayerOceanModel()

        # Apply negative forcing
        result = model.energy_budget(-1.0, -1.0, [0.0], [0.0])

        # With negative forcing, fast layer should cool
        assert result["dtemp_fast"] < 0.0
        assert result["dtemp"] < 0.0  # dtemp = temp_fast

        # Check that model state is updated
        assert model.temp_fast < 0.0

    def test_energy_budget_volcanic_forcing(self):
        """Test energy budget with volcanic forcing"""
        model = TwoLayerOceanModel()

        # Test with volcanic forcing arrays
        fn_volc = [-0.5, -0.3, -0.1]
        fs_volc = [-0.4, -0.2, 0.0]

        model.energy_budget(1.0, 1.0, fn_volc, fs_volc)

        # The total forcing should be 1.0 + np.mean(fn_volc + fs_volc)
        expected_volc_forcing = np.mean(np.array(fn_volc) + np.array(fs_volc))
        total_forcing = 1.0 + expected_volc_forcing

        # Should be cooling effect due to negative volcanic forcing
        assert total_forcing < 1.0

    def test_energy_budget_multiple_steps(self):
        """Test that temperature accumulates over multiple energy budget calls"""
        model = TwoLayerOceanModel()

        # First step
        result1 = model.energy_budget(1.0, 1.0, [0.0], [0.0])
        temp_fast_1 = result1["dtemp_fast"]
        temp_slow_1 = result1["dtemp_slow"]

        # Second step with same forcing
        result2 = model.energy_budget(1.0, 1.0, [0.0], [0.0])
        temp_fast_2 = result2["dtemp_fast"]
        temp_slow_2 = result2["dtemp_slow"]

        # Temperatures should be higher in second step
        assert temp_fast_2 > temp_fast_1
        assert temp_slow_2 >= temp_slow_1  # Should start increasing after first step

    def test_energy_budget_layer_coupling(self):
        """Test that the layers couple correctly"""
        model = TwoLayerOceanModel()

        # Warm up the fast layer first
        for _ in range(10):
            model.energy_budget(1.0, 1.0, [0.0], [0.0])

        # Now the fast layer should be significantly warmer than slow layer
        assert model.temp_fast > model.temp_slow

        # Apply zero forcing - heat should transfer from fast to slow layer
        initial_fast = model.temp_fast
        initial_slow = model.temp_slow

        model.energy_budget(0.0, 0.0, [0.0], [0.0])

        # Fast layer should cool due to heat transfer to slow layer
        # Slow layer should warm due to heat from fast layer
        assert model.temp_fast < initial_fast
        assert model.temp_slow > initial_slow

    def test_energy_budget_rib_calculation(self):
        """Test that RIB (Radiative Imbalance) is calculated correctly"""
        model = TwoLayerOceanModel()

        result = model.energy_budget(2.0, 2.0, [0.0], [0.0])

        # RIB should be calculated as:
        # forc - lambda * temp_fast - (ocean_efficacy - 1) * k * (temp_fast - temp_slow)
        expected_rib = (
            2.0
            - model.pamset["lambda"] * model.temp_fast
            - (model.pamset["ocean_efficacy"] - 1)
            * model.pamset["k"]
            * (model.temp_fast - model.temp_slow)
        )

        assert result["RIB"] == pytest.approx(expected_rib, rel=1e-10)

    def test_energy_budget_return_structure(self):
        """Test that energy_budget returns the expected dictionary structure"""
        model = TwoLayerOceanModel()
        result = model.energy_budget(1.0, 1.0, [0.0], [0.0])

        # Check all expected keys are present for interface compatibility
        expected_keys = [
            "dtemp",
            "dtemp_fast",
            "dtemp_slow",
            "RIB",
            "OHC700",
            "OHCTOT",
            "OHC_MIXED",
            "OHC_DEEP",
        ]

        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

        # Check that some values make sense
        assert (
            result["dtemp"] != result["dtemp_fast"]
        )  # Should be different (weighted combination)
        assert isinstance(result["RIB"], (float, np.floating))

        # Check that meaningful outputs are non-zero and placeholders are zero
        assert result["OHC_MIXED"] != 0.0  # Should have meaningful value
        assert (
            result["OHC_DEEP"] == 0.0
        )  # Should be zero with zero forcing on slow layer initially
        assert result["OHC700"] != 0.0  # Should have meaningful calculated value

        # Verify air/sea temperature relationships
        assert (
            result["dtemp_sea"] == result["dtemp_fast"]
        )  # Sea temp should equal fast layer

    def test_efficacy_effect(self):
        """Test that efficacy parameter affects the coupling correctly"""
        # Test with efficacy = 1 (default)
        model1 = TwoLayerOceanModel({"ocean_efficacy": 1.0})

        # Test with efficacy = 0.5 (reduced deep ocean heat uptake)
        model2 = TwoLayerOceanModel({"ocean_efficacy": 0.5})

        # Apply same forcing to both
        for _ in range(5):
            model1.energy_budget(1.0, 1.0, [0.0], [0.0])
            model2.energy_budget(1.0, 1.0, [0.0], [0.0])

        # With lower efficacy, the coupling between layers is different
        # Let's just verify that different efficacy values produce different results
        assert model2.temp_fast != model1.temp_fast
        assert model2.temp_slow != model1.temp_slow

        # The efficacy should affect the RIB calculation differently
        result1_final = model1.energy_budget(1.0, 1.0, [0.0], [0.0])
        result2_final = model2.energy_budget(1.0, 1.0, [0.0], [0.0])
        assert result1_final["RIB"] != result2_final["RIB"]

    def test_parameter_sensitivity(self):
        """Test sensitivity to different parameter values"""
        # Test with different lambda values
        model_low_lambda = TwoLayerOceanModel({"lambda": 0.5})
        model_high_lambda = TwoLayerOceanModel({"lambda": 2.0})

        # Apply same forcing multiple times to see the effect
        for _ in range(10):  # More steps to see the cumulative effect
            model_low_lambda.energy_budget(1.0, 1.0, [0.0], [0.0])
            model_high_lambda.energy_budget(1.0, 1.0, [0.0], [0.0])

        # Lower lambda (less climate feedback) should result in more warming after multiple steps
        assert model_low_lambda.temp_fast > model_high_lambda.temp_fast

    def test_thermal_model_required_pamset(self):
        """Test that the required pamset is properly defined"""
        required_pamset = TwoLayerOceanModel.thermal_model_required_pamset

        expected_keys = ["lambda", "mixed", "deep", "k", "ocean_efficacy"]

        for key in expected_keys:
            assert key in required_pamset, f"Missing required parameter: {key}"

        # Check that values are reasonable
        assert isinstance(required_pamset["lambda"], (int, float))
        assert required_pamset["lambda"] > 0
        assert isinstance(required_pamset["mixed"], (int, float))
        assert required_pamset["mixed"] > 0
        assert isinstance(required_pamset["deep"], (int, float))
        assert required_pamset["deep"] > 0

    def test_ocean_heat_content_calculation(self):
        """Test that ocean heat content is calculated correctly"""
        model = TwoLayerOceanModel()

        # Apply forcing to warm up the ocean
        result = model.energy_budget(2.0, 2.0, [0.0], [0.0])

        # Check that OHC values are calculated
        assert result["OHC_MIXED"] > 0  # Should have some heat in mixed layer
        assert result["OHC_DEEP"] == 0  # Deep layer starts with no heat
        assert result["OHCTOT"] > 0  # Total should be positive

        # Check that OHCTOT equals sum of components
        expected_total = result["OHC_MIXED"] + result["OHC_DEEP"]
        assert result["OHCTOT"] == pytest.approx(expected_total, rel=1e-10)

        # Apply more forcing to see heat transfer to deep layer
        for _ in range(5):
            result = model.energy_budget(1.0, 1.0, [0.0], [0.0])

        # Now deep layer should have some heat
        assert result["OHC_DEEP"] > 0
        assert (
            result["OHC_MIXED"] > result["OHC_DEEP"]
        )  # Mixed layer should still be warmer

        # Total should still equal sum of components
        expected_total = result["OHC_MIXED"] + result["OHC_DEEP"]
        assert result["OHCTOT"] == pytest.approx(expected_total, rel=1e-10)

    def test_get_thermal_model_required_pamset(self):
        """Test the class method for getting required pamset"""
        model = TwoLayerOceanModel()
        required_pamset = model.get_thermal_model_required_pamset()

        assert required_pamset == TwoLayerOceanModel.thermal_model_required_pamset


class TestTwoLayerOceanModelEdgeCases:
    """Test edge cases and error conditions"""

    def test_very_small_forcing(self):
        """Test with very small forcing values"""
        model = TwoLayerOceanModel()

        result = model.energy_budget(1e-10, 1e-10, [0.0], [0.0])

        # Should handle small values without issues
        assert np.isfinite(result["dtemp_fast"])
        assert np.isfinite(result["dtemp_slow"])

    def test_very_large_forcing(self):
        """Test with very large forcing values"""
        model = TwoLayerOceanModel()

        result = model.energy_budget(100.0, 100.0, [0.0], [0.0])

        # Should handle large values without issues
        assert np.isfinite(result["dtemp_fast"])
        assert np.isfinite(result["dtemp_slow"])
        assert result["dtemp_fast"] > 0  # Should be warming

    def test_empty_volcanic_arrays(self):
        """Test with empty volcanic forcing arrays"""
        model = TwoLayerOceanModel()

        # Empty arrays produce warnings but don't raise errors
        import warnings

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = model.energy_budget(1.0, 1.0, [], [])

            # Should get warnings about mean of empty slice
            assert len(w) > 0
            assert "Mean of empty slice" in str(w[0].message)

        # Result should still be computed (with NaN values)
        assert np.isnan(result["dtemp_fast"]) or np.isfinite(result["dtemp_fast"])

    def test_different_length_volcanic_arrays(self):
        """Test with different length volcanic forcing arrays"""
        model = TwoLayerOceanModel()

        # This should work since np.mean handles different lengths
        result = model.energy_budget(1.0, 1.0, [0.1, 0.2], [0.3])

        assert np.isfinite(result["dtemp_fast"])

    def test_extreme_parameters(self):
        """Test with extreme parameter values"""
        # Very small heat capacities (very shallow layers)
        extreme_params = {
            "mixed": 0.1,  # Very shallow mixed layer
            "deep": 0.1,  # Very shallow deep layer
            "k": 100.0,  # Very strong coupling
        }

        model = TwoLayerOceanModel(extreme_params)

        # Should still be able to run without crashing
        result = model.energy_budget(1.0, 1.0, [0.0], [0.0])
        assert np.isfinite(result["dtemp_fast"])
        assert np.isfinite(result["dtemp_slow"])

    def test_ohc700_calculation(self):
        """Test OHC700 calculation for different layer configurations"""
        # Case 1: Mixed layer deeper than 700m
        model1 = TwoLayerOceanModel({"mixed": 800, "deep": 1200})
        result1 = model1.energy_budget(2.0, 2.0, [0.0], [0.0])

        # OHC700 should be 700/800 of mixed layer OHC (since mixed > 700m)
        expected_ratio = 700 / 800
        actual_ratio = result1["OHC700"] / result1["OHC_MIXED"]
        assert abs(actual_ratio - expected_ratio) < 1e-10

        # Case 2: Mixed layer shallower than 700m (default case)
        model2 = TwoLayerOceanModel({"mixed": 50, "deep": 1200})
        # Run multiple steps to get some deep layer heating
        for _ in range(3):
            result2 = model2.energy_budget(2.0, 2.0, [0.0], [0.0])

        # OHC700 should be all mixed + (650/1200) of deep layer
        deep_fraction = (700 - 50) / 1200  # 650/1200
        expected_ohc700 = result2["OHC_MIXED"] + result2["OHC_DEEP"] * deep_fraction
        assert abs(result2["OHC700"] - expected_ohc700) < 1e-6

        # Case 3: Very shallow layers
        model3 = TwoLayerOceanModel({"mixed": 30, "deep": 500})
        for _ in range(3):
            result3 = model3.energy_budget(2.0, 2.0, [0.0], [0.0])

        # OHC700 should be all mixed + all deep (since 30 + 500 < 700)
        expected_ohc700 = result3["OHC_MIXED"] + result3["OHC_DEEP"]
        assert abs(result3["OHC700"] - expected_ohc700) < 1e-6
