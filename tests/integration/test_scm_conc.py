import os
import shutil

import numpy as np
import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM


def check_output(
    output_dir,
    expected_output_dir,
    update_expected_files=False,
    rtol=1e-2,
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

            pdt.assert_frame_equal(
                res.T,
                exp.T,
                check_like=True,
                rtol=rtol,
            )


def check_output_just_some_lines(
    output_dir,
    expected_output_dir,
    update_expected_files=False,
    rtol=1e-2,
    files=["output_temp.txt", "output_forc.txt"],
    lines=10,
):
    for filename in files:
        file_to_check = os.path.join(output_dir, filename)
        file_expected = os.path.join(expected_output_dir, filename)

        if update_expected_files:
            shutil.copyfile(file_to_check, file_expected)
        else:

            res = pd.read_csv(
                file_to_check, delim_whitespace=True, skiprows=range(lines, 352)
            )
            exp = pd.read_csv(
                file_expected, delim_whitespace=True, skiprows=range(lines, 352)
            )
            pdt.assert_index_equal(res.index, exp.index)

            pdt.assert_frame_equal(
                res.T,
                exp.T,
                check_like=True,
                rtol=rtol,
            )


def test_ciceroscm_run_emi(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run({"output_folder": outdir})

    check_output(outdir, os.path.join(test_data_dir, "ssp245_emis"))
    check_output_just_some_lines(
        outdir,
        os.path.join(test_data_dir, "ssp245_emis"),
        files=["output_forc.txt"],
        lines=19,
    )
    check_output_just_some_lines(
        outdir,
        os.path.join(test_data_dir, "ssp245_emis"),
        files=["output_temp.txt"],
        lines=16,
    )


def test_ciceroscm_short_run(tmpdir, test_data_dir):
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:
    nystart = 1900
    nyend = 2050
    emstart = 1950
    cscm = CICEROSCM(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nystart": nystart,
            "emstart": emstart,
            "nyend": nyend,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )

    cscm._run({"output_folder": outdir})

    file_results = os.path.join(outdir, "output_em.txt")
    exp_index = np.arange(nystart, nyend + 1)
    res = pd.read_csv(file_results, delim_whitespace=True)
    np.testing.assert_equal(res.Year.to_numpy(), exp_index)

    cscm._run({"results_as_dict": True})
    expected_keys = [
        "emissions",
        "concentrations",
        "forcing",
        "OHC700",
        "OHCTOT",
        "RIB_glob",
        "RIB_N",
        "RIB_S",
        "dT_glob",
        "dT_NH",
        "dT_SH",
        "dT_glob_air",
        "dT_NH_air",
        "dT_SH_air",
        "dT_glob_sea",
        "dT_NH_sea",
        "dT_SHsea",
        "dSL(m)",
        "dSL_thermal(m)",
        "dSL_ice(m)",
    ]
    for key in expected_keys:
        assert key in cscm.results
    # Put this in again, find out what is happening with CF4
    # check_output(
    #    outdir_save,
    #    os.path.join(test_data_dir, "ssp245_emis"),
    #    files=["output_temp.txt"],
    #    rtol=1.0,
    # )

    # check_output(
    #    outdir_save,
    #    os.path.join(test_data_dir, "ssp245_emis"),
    #    files=["output_forc.txt"],
    #    rtol=0.1,
    # )


def test_ciceroscm_run_conc(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "conc_run": True,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run({"output_folder": outdir})

    check_output(
        outdir,
        os.path.join(test_data_dir, "ssp245_conc"),
        files=["output_conc.txt", "output_em.txt", "output_forc.txt", "output_ohc.txt"],
    )

    """
    # 1pct CO2 without sunvolc

    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir,
        },
        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    )

    check_output(outdir, os.path.join(test_data_dir, "1pct_CO2_no_sunvolc"))

    #1 ppct CO2 with sunvolc
    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir, "sunvolc": 1,"nyend": 2100,
        },
        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    )

    check_output(outdir, os.path.join(test_data_dir, "1pct_CO2"))
    # check_output(outdir, os.path.join(test_data_dir,"1pct_CO2_no_sunvolc"))
    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir_save, "sunvolc": 1,"nyend": 2100,
            "threstemp": 0,
        },
        {"forc_file": os.path.join(test_data_dir, "CO2_1pros.txt")},
    )

    check_output_subset(outdir_save, os.path.join(test_data_dir, "nr_test_1pct_CO2"))
    #Test NR-setup:
    """


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
