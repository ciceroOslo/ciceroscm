import os

import numpy as np
import pandas as pd
import pytest

from ciceroscm import input_handler


def test_check_inputfiles():
    cfg = {"gaspamfile": "definitely_not_a_valid_path"}
    with pytest.raises(FileNotFoundError):
        cfg = input_handler.check_inputfiles(cfg)
    cfg = {"nat_ch4_file": "definitely_not_a_valid_path"}
    cfg = input_handler.check_inputfiles(cfg)
    assert cfg["nat_ch4_file"] == os.path.join(
        os.getcwd(), "input_OTHER", "NATEMIS", "natemis_ch4.txt"
    )
    cfg = {"concentration_file": "also_not_a_valid_path"}
    with pytest.raises(FileNotFoundError):
        cfg = input_handler.check_inputfiles(cfg)


def test_read_components(test_data_dir):
    df_gas = input_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    print(df_gas.columns)
    assert df_gas.columns.tolist() == [
        "EM_UNIT",
        "CONC_UNIT",
        "BETA",
        "ALPHA",
        "TAU1",
        "TAU2",
        "TAU3",
        "NAT_EM",
        "SARF_TO_ERF",
    ]
    assert len(df_gas.index) == 46


def test_read_forc(test_data_dir):
    forc = input_handler.read_forc(os.path.join(test_data_dir, "CO2_1pros.txt"))
    assert isinstance(forc, np.ndarray)
    assert len(forc) == 736

    forc = input_handler.read_forc(os.path.join(test_data_dir, "AR6_ERF_1750-2019.csv"))
    assert isinstance(forc, pd.DataFrame)
    total = forc["total"].to_numpy()
    assert len(total) == 270
    assert total[1] == 0.2866421077706115


def test_get_data_and_more(test_data_dir):
    ih = input_handler.InputHandler(
        {
            "forc_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
            "unicorn_emissions_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
        }
    )
    forc = ih.get_data("forc")
    assert isinstance(forc, np.ndarray)
    assert len(forc) == 736
    assert ih.optional_pam("forc")
    assert not ih.optional_pam("made_up_pam")
    assert not ih.conc_run()
    with pytest.raises(KeyError, match="No user or default data for emissions"):
        ih.get_data("emissions")
    with pytest.raises(KeyError, match="No user or default data for unicorn_emissions"):
        ih.get_data("unicorn_emissions")


def test_read_volc_sun(test_data_dir):
    ih = input_handler.InputHandler(
        {
            "forc_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
            "unicorn_emissions_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
        }
    )

    # Testing volcano rf reading
    df_volc = ih.read_data_on_year_row(
        os.path.join(test_data_dir, "meanVOLCmnd_ipcc_NH.txt")
    )
    print(df_volc.columns.tolist)
    assert len(df_volc.index) == 351
    assert df_volc.columns.tolist() == list(range(12))

    # Testing solar rf reading:
    df_sun = ih.read_data_on_year_row(os.path.join(test_data_dir, "solar_IPCC.txt"))
    assert len(df_sun.index) == 351
    assert df_sun.columns.tolist() == [0]

    # Testing LUCalbedo rf reading
    df_luc = ih.read_data_on_year_row(os.path.join(test_data_dir, "IPCC_LUCalbedo.txt"))
    assert len(df_luc.index) == 351
    assert df_luc.columns.tolist() == [0]
