import os

import numpy as np

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
            "fnpp_temp_coeff": 0,
        },
    )
    print("Attempting to change parameters in carbon_cycle_model")
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
            "fnpp_temp_coeff": 0,
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
            "fnpp_temp_coeff": 0,
        },
    )
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_emiconc={"fnpp_temp_coeff": -1, "mixed_carbon": 75.0},
    )
    decreased_fnpp_results = {
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
            "fnpp_temp_coeff": -1,
        },
    )
    print(decreased_fnpp_results["Biosphere carbon flux"])
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_emiconc={"fnpp_temp_coeff": 20},
    )
    increased_fnpp_results = {
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
            "fnpp_temp_coeff": 20,
        },
    )
    print(increased_fnpp_results["Biosphere carbon flux"])

    assert np.all(default_results["Temperature"] >= beta_f_results["Temperature"])
    assert np.all(
        default_results["Airborne fraction"] > beta_f_results["Airborne fraction"]
    )
    assert np.all(
        default_results["Biosphere carbon flux"]
        < beta_f_results["Biosphere carbon flux"]
    )
    assert np.all(
        default_results["CO2 concentration"] > beta_f_results["CO2 concentration"]
    )
    assert np.all(
        default_results["Ocean carbon flux"] > beta_f_results["Ocean carbon flux"]
    )

    for key, value in mixed_carbon_results.items():
        assert not np.isclose(value, default_results[key]).all()

    assert np.all(
        increased_fnpp_results["Temperature"] < decreased_fnpp_results["Temperature"]
    )
    assert np.all(
        default_results["Temperature"] < decreased_fnpp_results["Temperature"]
    )
    assert np.all(
        increased_fnpp_results["Biosphere carbon flux"]
        >= decreased_fnpp_results["Biosphere carbon flux"]
    )
    assert np.all(
        default_results["Biosphere carbon flux"]
        > decreased_fnpp_results["Biosphere carbon flux"]
    )
    assert np.all(
        increased_fnpp_results["CO2 concentration"]
        < decreased_fnpp_results["CO2 concentration"]
    )
    assert np.all(
        default_results["CO2 concentration"]
        < decreased_fnpp_results["CO2 concentration"]
    )
    assert np.all(
        increased_fnpp_results["Airborne fraction"]
        < decreased_fnpp_results["Airborne fraction"]
    )
    assert np.all(
        default_results["Airborne fraction"]
        < decreased_fnpp_results["Airborne fraction"]
    )
    assert np.all(
        increased_fnpp_results["Ocean carbon flux"]
        < decreased_fnpp_results["Ocean carbon flux"]
    )
    assert np.all(
        default_results["Ocean carbon flux"]
        < decreased_fnpp_results["Ocean carbon flux"]
    )
