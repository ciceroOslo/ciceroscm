import os

import numpy as np
import pandas as pd

from ciceroscm import CICEROSCM, input_handler, pub_utils


def check_output_equal_and_different_column(
    expected_output_dir,
    new_results,
    different=["Total_forcing"],
    equal=["CFC-11", "SO4_IND"],
    filename="output_forc.txt",
):
    exp = pd.read_csv(os.path.join(expected_output_dir, filename), sep="\t")
    print(exp.columns)
    for column in different:
        assert not np.allclose(exp[column], new_results[column])
    for column in equal:
        print(exp[column])
        print(new_results[column])
        assert np.allclose(exp[column], new_results[column], rtol=2e-3)


def test_ciceroscm_run_w_regional_aerosols(test_data_dir):
    df_gas = input_handler.read_components(
        os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt")
    )
    reg_aerosol_RF_file = os.path.join(test_data_dir, "HTAP_reg_aerosol_RF.txt")
    reg_aerosol_df = pd.read_csv(reg_aerosol_RF_file, sep="\t", index_col=0)
    reg_aerosol_df.rename(columns={"sulfate": "SO2"}, inplace=True)
    df_gas_updated = pub_utils.make_regional_aerosol_gaspamdata(df_gas, reg_aerosol_df)

    cscm = CICEROSCM(
        {
            # "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "gaspam_data": df_gas_updated,
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(
                test_data_dir, "ssp245_with_regional_aerosols_em_RCMIP.txt"
            ),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    assert cscm.ce_handler.pamset["qbc"] == 0
    assert cscm.ce_handler.pamset["qoc"] == 0
    assert cscm.ce_handler.pamset["qdirso2"] == 0
    cscm._run(
        {"results_as_dict": True},
        pamset_udm={
            "rlamdo": 15.1,
            "akapa": 0.657,
            "cpi": 0.208,
            "W": 2.2,
            "beto": 6.9,
            "lambda": 0.606,
            "mixed": 107.0,
        },
        pamset_emiconc={
            "qbmb": 0.0,
            "qo3": 0.5,
            "qdirso2": -0.3701 / 57,
            "qindso2": -0.4163 / 57,
            "qbc": 0.163 / 7,
            "qoc": -0.084 / 16,
            "qh2o_ch4": 0.171,
        },
    )
    assert cscm.ce_handler.pamset["qbc"] == 0
    assert cscm.ce_handler.pamset["qoc"] == 0
    assert cscm.ce_handler.pamset["qdirso2"] == 0
    for region in ["ASIA", "LAM", "MAF", "OECD", "REF"]:
        for aerosol in ["SO2", "BC", "OC"]:
            assert f"{aerosol}_{region}" in cscm.results["emissions"].keys()
            assert f"{aerosol}_{region}" in cscm.results["forcing"].keys()
    print(cscm.results["forcing"].keys())
    assert np.allclose(cscm.results["forcing"]["SO4_DIR"], 0)
    assert np.allclose(cscm.results["forcing"]["BC"], 0)
    assert np.allclose(cscm.results["forcing"]["OC"], 0)
    check_output_equal_and_different_column(
        os.path.join(test_data_dir, "ssp245_emis"), cscm.results["forcing"]
    )
    # TODO: Check output differs from standard, but only for aerosols
