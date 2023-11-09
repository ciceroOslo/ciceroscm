# Tests for different decay functions will go here...
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM, pub_utils


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
            res = pd.read_csv(file_to_check, delim_whitespace=True)
            exp = pd.read_csv(file_expected, delim_whitespace=True)
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
                "rs_function": lambda it, idtm: 0.5 + 0.5 * np.exp(-it / idtm / 100.0),
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
    rb_decay = pub_utils.make_rb_function_from_arrays(
        [0.70211, 13.4141e-3, -0.71846, 2.9323e-3], [1 / 0.35, 20.0, 120 / 55.0, 100.0]
    )
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
            "rb_function": rb_decay,
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
    rb_decay = pub_utils.make_rb_function_from_arrays(
        [0.70211, 13.4141e-3, -0.71846, 2.9323e-3], [1 / 0.35, 20.0, 120 / 55.0, 100.0]
    )
    rs_decay_fair = pub_utils.make_rs_function_from_arrays(
        [0.24278, 0.13963, 0.089318, 0.03782, 0.035549],
        [1.2679, 5.2528, 18.601, 68.736, 232.3],
    )
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
            "rb_function": rb_decay,
            "rs_function": rs_decay_fair,
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
