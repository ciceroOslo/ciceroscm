import os
import re

import numpy as np
import pandas as pd
import pytest

from ciceroscm import input_handler, pub_utils
from ciceroscm.carbon_cycle import rfuns


def test_making_biotic_decay_function():
    rb_C = np.array([0.5, 0.25, 0.25])
    rb_T = np.array([2.5, 10, 60])
    rb_C, rb_T = pub_utils._check_array_consistency(rb_C, rb_T)
    rb_func = rfuns.rb_function2(np.arange(25), rb_coef=rb_C, rb_tim=rb_T, idtm=24)
    expected = (
        0.70211 * np.exp(-0.35)
        + 13.4141e-3 * np.exp(-1 / 20.0)
        - 0.71846 * np.exp(-55.0 / 120)
        + 2.9323e-3 * np.exp(-1 / 100.0)
    )
    assert rb_func[0] == 0.0
    assert np.allclose(rb_func[24], expected, rtol=1e-2)
    # assert (rb_func(24, 24) - expected) / expected < 1.0e-9


def test_making_carbon_pool_decay_function():
    rs_C = np.array([0.24278, 0.13963, 0.089318, 0.03782, 0.035549])
    rs_T = np.array([5.2528, 18.601, 68.736, 232.3])
    rs_C, rs_T = pub_utils._check_array_consistency(rs_C, rs_T, for_rs=True)
    rs_func = rfuns.rs_function2(np.arange(1001), rs_coef=rs_C, rs_tim=rs_T, idtm=24)
    expected = (
        0.24278
        + 0.13963 * np.exp(-1 / 5.2528)
        + 0.089318 * np.exp(-1 / 18.601)
        + 0.03782 * np.exp(-1 / 68.736)
        + 0.035549 * np.exp(-1 / 232.3)
    ) / np.sum(rs_C)
    assert rs_func[0] == 1.0
    assert rs_func[1000] < 1.0
    assert (rs_func[24] - expected) / expected < 1e-7


def test_error_management_decay_functions():
    rs_T = "Not an array"
    rs_C = [1, 4, -5]
    with pytest.raises(TypeError):
        pub_utils._check_array_consistency(rs_C, rs_T)
    rs_T = [3, 4]
    with pytest.raises(
        ValueError,
        match="Coefficient array length must be equal to exponent array length",
    ):
        pub_utils._check_array_consistency(rs_C, rs_T)
    with pytest.raises(ValueError, match="All timescales must be positive numbers"):
        pub_utils._check_array_consistency(rs_C, rs_C)
    with pytest.raises(ValueError, match="All timescales must be positive numbers"):
        pub_utils._check_array_consistency(rs_C, rs_C)
    rs_T = [1, 20, 100]
    with pytest.raises(ValueError, match="All coefficients must be positive"):
        pub_utils._check_array_consistency(rs_C, rs_T)
    with pytest.raises(
        ValueError,
        match=re.escape(
            "Coefficient array length must be equal to exponent array length + 1"
        ),
    ):
        pub_utils._check_array_consistency(rs_T, rs_T, for_rs=True)


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


def test_make_regional_aerosol_gaspamdata(test_data_dir):
    df_gas = input_handler.read_components(
        os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt")
    )
    reg_aerosol_RF_file = os.path.join(test_data_dir, "HTAP_reg_aerosol_RF.txt")
    reg_aerosol_df = pd.read_csv(reg_aerosol_RF_file, sep="\t", index_col=0)
    reg_aerosol_df.rename(columns={"sulfate": "SO2"}, inplace=True)
    df_gas_updated = pub_utils.make_regional_aerosol_gaspamdata(df_gas, reg_aerosol_df)
    aerosols = reg_aerosol_df.columns.tolist()
    regions = pub_utils.REG_MAPPING_DEFAULT
    assert (
        df_gas.shape[0] + len(aerosols) * len(regions.keys()) == df_gas_updated.shape[0]
    )
    for region, reg_map in regions.items():
        for aerosol in aerosols:
            compound_name = f"{aerosol}_{region}"
            assert compound_name in df_gas_updated.index
            if compound_name in df_gas_updated.index:
                assert (
                    df_gas_updated.loc[compound_name, "ALPHA"]
                    == reg_aerosol_df.loc[reg_map, aerosol] * 1e9
                )
