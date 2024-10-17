import os

import numpy as np

# import matplotlib.pyplot as plt
from ciceroscm import CICEROSCM


def check_post_run_parameter_values(cscm, dict_values):
    for key, item in dict_values.items():
        assert cscm.ce_handler.carbon_cycle.pamset[key] == item


def test_changing_carbon_cycle_parameters(test_data_dir):
    nystart = 1900
    nyend = 2015
    emstart = 1950
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nystart": nystart,
            "emstart": emstart,
            "nyend": nyend,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "idtm": 24,
        },
    )
    cscm._run({"results_as_dict": True, "carbon_cycle_outputs": True})
    default_results = {
        "Biosphere carbon flux": cscm.results["carbon cycle"][
            "Biosphere carbon flux"
        ].values.copy(),
        "Ocean carbon flux": cscm.results["carbon cycle"][
            "Ocean carbon flux"
        ].values.copy(),
        "Airborne fraction": cscm.results["carbon cycle"][
            "Airborne fraction CO2"
        ].values.copy(),
        "Temperature": cscm.results["dT_glob"].copy(),
        "CO2 concentration": cscm.results["concentrations"]["CO2"].values.copy(),
    }
    check_post_run_parameter_values(
        cscm,
        {
            "beta_f": 0.287,
            "mixed_carbon": 75.0,
        },
    )
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_emiconc={"beta_f": 0.3},
    )
    beta_f_results = {
        "Biosphere carbon flux": cscm.results["carbon cycle"][
            "Biosphere carbon flux"
        ].values.copy(),
        "Ocean carbon flux": cscm.results["carbon cycle"][
            "Ocean carbon flux"
        ].values.copy(),
        "Airborne fraction": cscm.results["carbon cycle"][
            "Airborne fraction CO2"
        ].values.copy(),
        "Temperature": cscm.results["dT_glob"].copy(),
        "CO2 concentration": cscm.results["concentrations"]["CO2"].values.copy(),
    }
    check_post_run_parameter_values(
        cscm,
        {
            "beta_f": 0.3,
            "mixed_carbon": 75.0,
        },
    )
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_emiconc={"mixed_carbon": 107.0, "beta_f": 0.287},
    )
    mixed_carbon_results = {
        "Biosphere carbon flux": cscm.results["carbon cycle"][
            "Biosphere carbon flux"
        ].values.copy(),
        "Ocean carbon flux": cscm.results["carbon cycle"][
            "Ocean carbon flux"
        ].values.copy(),
        "Airborne fraction": cscm.results["carbon cycle"][
            "Airborne fraction CO2"
        ].values.copy(),
        "Temperature": cscm.results["dT_glob"].copy(),
        "CO2 concentration": cscm.results["concentrations"]["CO2"].values.copy(),
    }
    check_post_run_parameter_values(
        cscm,
        {
            "beta_f": 0.287,
            "mixed_carbon": 107.0,
        },
    )

    for key, value in beta_f_results.items():
        if key == "Biosphere carbon flux":
            assert np.all(default_results[key] <= value)
        else:
            assert np.all(default_results[key] >= value)

    for key, value in mixed_carbon_results.items():
        if key == "Ocean carbon flux":
            assert np.all(default_results[key] <= value)
        else:
            assert np.all(default_results[key] >= value)
