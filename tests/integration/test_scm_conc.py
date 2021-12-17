import os
import shutil

import pandas as pd
import pandas.testing as pdt

from ciceroscm import CICEROSCM


def check_output(
    output_dir, expected_output_dir, update_expected_files=False, rtol=1e-2
):
    files = ["output_em.txt", "output_conc.txt","output_forc.txt","output_temp.txt", "output_ohc.txt"]

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


def test_ciceroscm_run_emi(tmpdir, test_data_dir):
    cscm = CICEROSCM()
    outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    #One year forcing:
    
    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir_save,
            "nyend": 2100,
        },
        {"concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"), "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"), "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"), "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt")},
    )
    
    check_output(outdir_save, os.path.join(test_data_dir, "ssp245_emis"))
"""
def test_ciceroscm_run_conc(tmpdir, test_data_dir):
    cscm = CICEROSCM()
    outdir_save = os.path.join(os.getcwd(), "output")
    outdir = str(tmpdir)
    #One year forcing:
    
    cscm._run(
        {
            "gaspamfile": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "output_prefix": outdir_save,
            "nyend": 2100,
        },
        {"concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"), "conc_run":True, "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"), "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"), "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt")},
    )
    
    check_output(outdir, os.path.join(test_data_dir, "1_year_blipp"))
   

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
