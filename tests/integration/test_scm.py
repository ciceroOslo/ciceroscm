import os
import shutil

import pandas as pd
import pandas.testing as pdt


from ciceroscm import CICEROSCM


def check_output(
    output_dir, expected_output_dir, update_expected_files=False, rtol=1e-2
):
    files = ["output_temp.txt", "output_rib.txt", "output_temp.txt"]

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
                res.T, exp.T, check_like=True, rtol=rtol,
            )


"""
def test_ciceroscm_zero_run(test_data_dir):
    cscm = CICEROSCM()
    outdir = os.path.join(os.getcwd(), "output")
    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir,
        },
        {"forc_file": os.path.join(test_data_dir, "zero_forcing.txt")},
    )

    check_output(outdir, os.path.join(test_data_dir, "all_zero"))
    # check_output(outdir, os.path.join(test_data_dir,"1pct_CO2_no_sunvolc"))

"""


def test_ciceroscm_run(tmpdir, test_data_dir):
    cscm = CICEROSCM()
    # outdir = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir,
            "nyend": 2100,
        },
        {"forc_file": os.path.join(test_data_dir, "test_forcing.txt")},
    )

    check_output(outdir, os.path.join(test_data_dir, "1_year_blipp"))

    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir,
        },
        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    )

    check_output(outdir, os.path.join(test_data_dir, "1pct_CO2_no_sunvolc"))
    # check_output(outdir, os.path.join(test_data_dir,"1pct_CO2_no_sunvolc"))


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
