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


def test_gaspam_compatibility_check(test_data_dir):
    ih = input_handler.InputHandler(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
        }
    )
    gaspam_testing = ih.get_data("gaspam")
    gaspam_testing.loc["CH4", "EM_UNIT"] = "Pg_C"
    gaspam_testing.loc["N2O", "CONC_UNIT"] = "NOUNIT"
    with pytest.raises(ValueError, match=r"Emissions file *"):
        input_handler.InputHandler(
            {
                "gaspam_data": gaspam_testing,
                "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            }
        )
    with pytest.raises(ValueError, match=r"Concentrations file *"):
        input_handler.InputHandler(
            {
                "gaspam_data": gaspam_testing,
                "concentrations_file": os.path.join(
                    test_data_dir, "ssp245_conc_RCMIP.txt"
                ),
            }
        )


def test_read_forc_with_aerosol_pattern_effect(test_data_dir, tmpdir):
    # Case 1: CO2, BC, SO4_IND (SO4_IND has negative values)
    forc_csv_1 = tmpdir.join("forc_aero_basic.csv")
    forc_csv_1.write(
        "year,CO2,BC,SO4_IND\n"
        "2000,2.0,0.4,-0.8\n"
        "2001,3.0,0.3,-0.3\n"
        "2002,1.0,0.5,-0.5\n"
    )
    df1 = input_handler.read_forc(str(forc_csv_1))
    assert "w_aero" in df1.columns
    # w_aero = (|BC| + |SO4_IND|) / (|CO2| + |BC| + |SO4_IND|)
    expected_w_aero_1 = np.array(
        [
            (0.4 + 0.8) / (2.0 + 0.4 + 0.8),  # 1.2 / 3.2 = 0.375
            (0.3 + 0.3) / (3.0 + 0.3 + 0.3),  # 0.6 / 3.6 = 1/6
            (0.5 + 0.5) / (1.0 + 0.5 + 0.5),  # 1.0 / 2.0 = 0.5
        ]
    )
    np.testing.assert_allclose(df1["w_aero"].to_numpy(), expected_w_aero_1)

    # Case 2: CO2, BC, SO4_IND, OC_Asia, OC_Africa (all aerosol cols negative except BC)
    # OC_Asia and OC_Africa start with "OC" so they are picked up by aero_cols_reg
    forc_csv_2 = tmpdir.join("forc_aero_regional_oc.csv")
    forc_csv_2.write(
        "year,CO2,BC,SO4_IND,OC_Asia,OC_Africa\n"
        "2000,2.0,0.4,-0.8,-0.2,-0.1\n"
        "2001,3.0,0.3,-0.3,-0.15,-0.05\n"
        "2002,1.0,0.5,-0.5,-0.1,-0.1\n"
    )
    df2 = input_handler.read_forc(str(forc_csv_2))
    assert "w_aero" in df2.columns
    # w_aero = (|BC| + |SO4_IND| + |OC_Asia| + |OC_Africa|) / (|CO2| + |BC| + |SO4_IND| + |OC_Asia| + |OC_Africa|)
    expected_w_aero_2 = np.array(
        [
            (0.4 + 0.8 + 0.2 + 0.1) / (2.0 + 0.4 + 0.8 + 0.2 + 0.1),  # 1.5 / 3.5
            (0.3 + 0.3 + 0.15 + 0.05) / (3.0 + 0.3 + 0.3 + 0.15 + 0.05),  # 0.8 / 3.8
            (0.5 + 0.5 + 0.1 + 0.1) / (1.0 + 0.5 + 0.5 + 0.1 + 0.1),  # 1.2 / 2.2
        ]
    )
    np.testing.assert_allclose(df2["w_aero"].to_numpy(), expected_w_aero_2)

    # Case 3: same columns as Case 1 but w_aero is already present in the file;
    # read_forc should pass through the precalculated values unchanged.
    precalc_w_aero = expected_w_aero_1
    forc_csv_3 = tmpdir.join("forc_aero_precalc_w_aero.csv")
    forc_csv_3.write(
        f"year,CO2,BC,SO4_IND,w_aero\n"
        f"2000,2.0,0.4,-0.8,{precalc_w_aero[0]}\n"
        f"2001,3.0,0.3,-0.3,{precalc_w_aero[1]}\n"
        f"2002,1.0,0.5,-0.5,{precalc_w_aero[2]}\n"
    )
    df3 = input_handler.read_forc(str(forc_csv_3))
    assert "w_aero" in df3.columns
    np.testing.assert_allclose(df3["w_aero"].to_numpy(), precalc_w_aero)
