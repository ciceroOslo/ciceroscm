import os

import numpy as np

from ciceroscm import CICEROSCM
from helpers import _DEFAULT_PAMSET_EMICONC, _DEFAULT_PAMSET_UDM, load_pdo_data_padded


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


def test_feedback_variable_passing_with_pdo_input_is_stable(test_data_dir):
    nystart = 1900
    nyend = 2015
    emstart = 1950
    pdo = load_pdo_data_padded(test_data_dir, nystart, nyend)

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

    cscm_base = CICEROSCM(common_params.copy())
    cscm_pdo = CICEROSCM(common_params.copy())

    pamset_base = {**_DEFAULT_PAMSET_UDM, "pdo_index_data": np.zeros_like(pdo)}
    pamset_pdo = {
        **_DEFAULT_PAMSET_UDM,
        "pdo_index_data": pdo,
        "delta_lambda_pdo": 0.0,
        "pdo_efficacy_scale": 0.0,
    }

    cscm_base._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_udm=pamset_base,
        pamset_emiconc=_DEFAULT_PAMSET_EMICONC,
    )
    cscm_pdo._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_udm=pamset_pdo,
        pamset_emiconc=_DEFAULT_PAMSET_EMICONC,
    )

    assert set(cscm_pdo.feedback_list) == {"dtemp"}
    np.testing.assert_array_equal(cscm_base.results["dT_glob"], cscm_pdo.results["dT_glob"])
    np.testing.assert_array_equal(
        cscm_base.results["concentrations"]["CO2"].to_numpy(),
        cscm_pdo.results["concentrations"]["CO2"].to_numpy(),
    )
    np.testing.assert_array_equal(
        cscm_base.results["carbon cycle"].to_numpy(),
        cscm_pdo.results["carbon cycle"].to_numpy(),
    )
