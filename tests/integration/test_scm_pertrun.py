import os
import shutil

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


def test_ciceroscm_run_pert_forc(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "perturb_forc_file": os.path.join(test_data_dir, "pertforc_test.txt"),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_emiconc={
            "qh2o_ch4": 0.171,
            "qbmb": 0.03,
            "qo3": 0.4,
            "qdirso2": -0.457,
            "qindso2": -0.514,
            "qbc": 0.200,
            "qoc": -0.103,
            "qh2o_ch4": 0.171,
        },
    )

    check_output(
        outdir, os.path.join(test_data_dir, "pert_tests"), files=["output_forc.txt"]
    )
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "perturb_forc_data": pd.read_csv(
                os.path.join(test_data_dir, "pertforc_test.txt")
            ),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run(
        {"output_folder": outdir},
        pamset_emiconc={
            "qh2o_ch4": 0.171,
            "qbmb": 0.03,
            "qo3": 0.4,
            "qdirso2": -0.457,
            "qindso2": -0.514,
            "qbc": 0.200,
            "qoc": -0.103,
            "qh2o_ch4": 0.171,
        },
    )

    check_output(
        outdir, os.path.join(test_data_dir, "pert_tests"), files=["output_forc.txt"]
    )


def test_ciceroscm_run_pert_emis(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "perturb_em_file": os.path.join(test_data_dir, "pertem_test.txt"),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run({"output_folder": outdir})

    check_output(
        outdir,
        os.path.join(test_data_dir, "pert_tests"),
        files=["output_em.txt"],
    )
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "perturb_em_data": pd.read_csv(
                os.path.join(test_data_dir, "pertem_test.txt"), index_col=None
            ),
        },
    )
    # outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    # One year forcing:

    cscm._run({"output_folder": outdir})

    check_output(
        outdir,
        os.path.join(test_data_dir, "pert_tests"),
        files=["output_em.txt"],
    )
