import os

import numpy as np
import pytest

from ciceroscm import input_handler, pub_utils


def test_making_biotic_decay_function():
    rb_C = np.array([0.70211, 13.4141e-3, -0.71846, 2.9323e-3])
    rb_T = np.array([1 / 0.35, 20.0, 120 / 55.0, 100.0])
    rb_func = pub_utils.make_rb_function_from_arrays(rb_C, rb_T)
    expected = (
        0.70211 * np.exp(-0.35)
        + 13.4141e-3 * np.exp(-1 / 20.0)
        - 0.71846 * np.exp(-55.0 / 120)
        + 2.9323e-3 * np.exp(-1 / 100.0)
    )
    print(rb_func(24, 24))
    print(expected)
    assert rb_func(0, 24) == np.sum(rb_C)
    assert rb_func(24, 24) == np.dot(rb_C, np.exp(-1.0 / rb_T))
    assert (rb_func(24, 24) - expected) / expected < 1.0e-9


def test_making_carbon_pool_decay_function():
    rs_C = np.array([0.24278, 0.13963, 0.089318, 0.03782, 0.035549])
    rs_T = np.array([1.2679, 5.2528, 18.601, 68.736, 232.3])
    rs_func = pub_utils.make_rs_function_from_arrays(rs_C, rs_T)
    expected = (
        (1 - np.sum(rs_C))
        + 0.24278 * np.exp(-1 / 1.2679)
        + 0.13963 * np.exp(-1 / 5.2528)
        + 0.089318 * np.exp(-1 / 18.601)
        + 0.03782 * np.exp(-1 / 68.736)
        + 0.035549 * np.exp(-1 / 232.3)
    )
    assert rs_func(0, 24) == 1.0
    assert rs_func(1000, 24) < 1.0
    print(expected)
    print(rs_func(24, 24))
    assert rs_func(24, 24) == (1 - np.sum(rs_C)) + np.dot(rs_C, np.exp(-1.0 / rs_T))
    assert (rs_func(24, 24) - expected) / expected < 1e-7


def test_error_management_decay_functions():
    rs_T = "Not an array"
    rs_C = [1, 4, -5]
    with pytest.raises(TypeError):
        pub_utils.make_rs_function_from_arrays(rs_C, rs_T)
    rs_T = [3, 4]
    with pytest.raises(
        ValueError,
        match="Coefficient array length must be equal to exponent array length",
    ):
        rs = pub_utils.make_rs_function_from_arrays(rs_C, rs_T)
    with pytest.raises(ValueError, match="All timescales must be positive numbers"):
        rb = pub_utils.make_rb_function_from_arrays(rs_C, rs_C)
    with pytest.raises(ValueError, match="All timescales must be positive numbers"):
        rs = pub_utils.make_rs_function_from_arrays(rs_C, rs_C)
    rs_T = [1, 20, 100]
    with pytest.raises(ValueError, match="All coefficients must be positive"):
        rs = pub_utils.make_rs_function_from_arrays(rs_C, rs_T)
    rb = pub_utils.make_rb_function_from_arrays(rs_C, rs_T)
    assert rb(0, 100) == np.sum(rs_C)
    with pytest.raises(
        ValueError,
        match="Sum of carbon decay coefficients must be less than or equal to 1",
    ):
        rs = pub_utils.make_rs_function_from_arrays(rs_T, rs_T)
    rs_C = [0.25, 0.25, 0.25]
    rs = pub_utils.make_rs_function_from_arrays(rs_C, rs_T)
    assert rs(0, 5364) == 1.0


def test_find_num_cl_in_hcfc():
    assert pub_utils.find_num_cl_in_hcfc("CFC-11") == 3
    assert pub_utils.find_num_cl_in_hcfc("CFC-12") == 2
    assert pub_utils.find_num_cl_in_hcfc("CFC-113") == 3
    assert pub_utils.find_num_cl_in_hcfc("CFC-114") == 2
    assert pub_utils.find_num_cl_in_hcfc("CFC-115") == 1
    assert pub_utils.find_num_cl_in_hcfc("HCFC-22") == 1
    assert pub_utils.find_num_cl_in_hcfc("HCFC-141b") == 2
    assert pub_utils.find_num_cl_in_hcfc("HCFC-123") == 2
    assert pub_utils.find_num_cl_in_hcfc("HCFC-142b") == 1


def test_make_cl_and_br_dictionaries(test_data_dir):
    df_gas = input_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    chlor_dict, brom_dict = pub_utils.make_cl_and_br_dictionaries(df_gas.index)
    assert brom_dict == {"H-1211": 1, "H-1301": 1, "CH3Br": 1, "H-2402": 2}
    chlor_dict_old = {
        "CFC-11": 3,
        "CFC-12": 2,
        "CFC-113": 3,
        "CFC-114": 2,
        "CFC-115": 1,
        "CCl4": 4,
        "CH3CCl3": 3,
        "HCFC-22": 1,
        "HCFC-141b": 2,
        "HCFC-123": 2,
        "HCFC-142b": 1,
    }
    assert chlor_dict != chlor_dict_old
    chlor_dict_old["H-1211"] = 1
    assert chlor_dict == chlor_dict_old
