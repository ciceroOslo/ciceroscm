import os
import shutil

import numpy as np
import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM, input_handler


def check_output(
    output_dir, expected_output_dir, update_expected_files=False, rtol=1e-2
):
    files = ["output_temp.txt", "output_ohc.txt", "output_sunvolc.txt"]

    for filename in files:
        file_to_check = os.path.join(output_dir, filename)
        file_expected = os.path.join(expected_output_dir, filename)

        if update_expected_files:
            shutil.copyfile(file_to_check, file_expected)

        else:
            res = pd.read_csv(file_to_check, delim_whitespace=True)
            exp = pd.read_csv(file_expected, delim_whitespace=True)
            pdt.assert_index_equal(res.index, exp.index)

            pdt.assert_frame_equal(
                res.T,
                exp.T,
                check_like=True,
                rtol=rtol,
            )


def check_output_subset(
    output_dir, expected_output_dir, update_expected_files=False, rtol=1e-2
):
    files = {"output_temp.txt": [0, 1, 2, 3], "output_ohc.txt": [0, 1, 2]}

    for filename in files:
        file_to_check = os.path.join(output_dir, filename)
        file_expected = os.path.join(expected_output_dir, filename)

        if update_expected_files:
            shutil.copyfile(file_to_check, file_expected)

        else:
            res = pd.read_csv(
                file_to_check, delim_whitespace=True, usecols=files[filename]
            )
            exp = pd.read_csv(file_expected, delim_whitespace=True)
            pdt.assert_index_equal(res.index, exp.index)

            pdt.assert_frame_equal(
                res.T,
                exp.T,
                check_like=True,
                rtol=rtol,
            )


def test_ciceroscm_run(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "forc_data": np.loadtxt(os.path.join(test_data_dir, "test_forcing.txt")),
        },
    )
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
            "ldtime": 12,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "1_year_blipp"))

    # 1pct CO2 without sunvolc
    cscm = CICEROSCM(
        {
            "nyend": 2100,
            "forc_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
        },
    )

    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "1pct_CO2_no_sunvolc"))

    ih = input_handler.InputHandler({"nyend": 2100, "sunvolc": 1})
    # 1 ppct CO2 with sunvolc
    cscm = CICEROSCM(
        {
            "nyend": 2100,
            "sunvolc": 1,
            "forc_file": os.path.join(test_data_dir, "CO2_1pros.txt"),
            "rf_volc_n_data": ih.get_data("rf_volc_n") + 0.371457071,
            "rf_volc_s_data": ih.get_data("rf_volc_s") + 0.353195076,
        },
    )
    outdir = str(tmpdir)
    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "1pct_CO2"))
    # check_output(outdir, os.path.join(test_data_dir,"1pct_CO2_no_sunvolc"))
    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "threstemp": 0,
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
        },
    )

    check_output_subset(outdir, os.path.join(test_data_dir, "nr_test_1pct_CO2"))
    # Test hemispheric split:
    cscm = CICEROSCM(
        {
            "nyend": 2100,
            "forc_file": os.path.join(test_data_dir, "test_forcing_hemisplit.txt"),
        },
    )
    outdir = str(tmpdir)
    # outdir = os.path.join(os.getcwd(), "output")
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "threstemp": 0,
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "1_year_blipp"))
    # Test component split:
    cscm = CICEROSCM(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "forc_file": os.path.join(test_data_dir, "test_forcing_components.txt"),
        },
    )
    outdir = str(tmpdir)
    # outdir = os.path.join(os.getcwd(), "output")
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_udm={
            "threstemp": 0,
            "rlamdo": 16.0,
            "akapa": 0.634,
            "cpi": 0.4,
            "W": 4.0,
            "beto": 3.5,
            "lambda": 0.540,
            "mixed": 60.0,
        },
    )

    check_output(outdir, os.path.join(test_data_dir, "1_year_blipp"))


"""
def test_cfg(test_data_dir):
    cscm = CICEROSCM()
    # cscm._run(
    #    {"gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),},
    #    {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    # )


#    cscm._run(
#        {"gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"), "sunvolc": 1},
#        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
#    )
"""
