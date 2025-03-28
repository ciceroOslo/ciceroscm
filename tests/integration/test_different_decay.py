# Tests for different decay functions will go here...
import os
import shutil
import warnings

import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM


def check_output_not_equal(
    output_dir,
    expected_output_dir,
    update_expected_files=False,
    rtol=1e-1,
    files=["output_em.txt", "output_conc.txt"],
):
    for filename in files:
        file_to_check = os.path.join(output_dir, filename)
        file_expected = os.path.join(expected_output_dir, filename)

        if update_expected_files:
            shutil.copyfile(file_to_check, file_expected)
        else:
            res = pd.read_csv(file_to_check, sep=r"\s+")
            exp = pd.read_csv(file_expected, sep=r"\s+")
            pdt.assert_index_equal(res.index, exp.index)

            try:
                pdt.assert_frame_equal(
                    res.T,
                    exp.T,
                    check_like=True,
                    rtol=rtol,
                )
            except AssertionError:
                pass
            else:
                raise AssertionError("DataFrames that should be different are equal")


def test_swap_out_rs_function(tmpdir, test_data_dir):
    with warnings.catch_warnings():
        cscm = CICEROSCM(
            {
                "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
                # "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
                "nyend": 2100,
                "nystart": 1750,
                "emstart": 1850,
                "concentrations_file": os.path.join(
                    test_data_dir, "ssp245_conc_RCMIP.txt"
                ),
                "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
                "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
                "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
                "rs_function": {
                    "coeffs": [0.5, 0.5],
                    "timescales": [100.0],
                },  # lambda it, idtm: 0.5 + 0.5 * np.exp(-it / idtm / 100.0),
            },
        )
        warnings.simplefilter("error")
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:
    cscm._run(
        {"output_folder": outdir},
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
            "qdirso2": -0.3701,
            "qindso2": -0.4163,
            "qbc": 0.163,
            "qoc": -0.084,
            "qh2o_ch4": 0.171,
        },
    )

    check_output_not_equal(
        outdir,
        os.path.join(test_data_dir, "ssp245_conc"),
        files=["output_conc.txt"],
    )


def test_swap_out_rb_function(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            # "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "rb_function": {"coeffs": [0.5, 0.25, 0.25], "timescales": [2.5, 10, 60]},
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:
    cscm._run(
        {"output_folder": outdir},
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
            "qdirso2": -0.3701,
            "qindso2": -0.4163,
            "qbc": 0.163,
            "qoc": -0.084,
            "qh2o_ch4": 0.171,
        },
    )
    check_output_not_equal(
        outdir,
        os.path.join(test_data_dir, "ssp245_conc"),
        files=["output_conc.txt"],
    )


def test_swap_out_both_r_functions(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            # "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "rb_function": {"coeffs": [0.5, 0.25, 0.25], "timescales": [2.5, 10, 60]},
            "rs_function": {
                "coeffs": [0.24278, 0.13963, 0.089318, 0.03782, 0.035549],
                "timescales": [5.2528, 18.601, 68.736, 232.3],
            },
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:
    cscm._run(
        {"output_folder": outdir},
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
            "qdirso2": -0.3701,
            "qindso2": -0.4163,
            "qbc": 0.163,
            "qoc": -0.084,
            "qh2o_ch4": 0.171,
        },
    )

    check_output_not_equal(
        outdir,
        os.path.join(test_data_dir, "ssp245_conc"),
        files=["output_conc.txt"],
    )
