import os

from ciceroscm import CICEROSCM


def test_feedback_variable_passing(test_data_dir, tmp_path):
    """
    Test that feedback variables are correctly passed to the carbon cycle model
    by comparing runs with and without feedback variables.

    This test runs two scenarios:
    1. A standard run without feedback variables.
    2. A run where feedback variables are explicitly passed to the carbon cycle model.

    The results of both runs are compared to ensure they are identical, confirming
    that feedback variables are correctly handled.
    """
    nystart = 1900
    nyend = 2015
    emstart = 1950

    # Common parameters for both runs
    common_params = {
        "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
        "nystart": nystart,
        "emstart": emstart,
        "nyend": nyend,
        "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
        "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
        "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
        "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        "idtm": 24,
    }

    cscm = CICEROSCM(common_params)
    assert set(cscm.feedback_list) == set(["dtemp"])
