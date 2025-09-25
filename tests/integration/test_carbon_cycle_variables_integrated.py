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
    print("Attempting to change parameters in carbon_cycle_model")

    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_carbon={"beta_f": 0.3},
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
        pamset_carbon={"mixed_carbon": 107.0, "beta_f": 0.287},
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
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_carbon={"solubility_sens": 0.03, "t_half": 0.2},
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
            "mixed_carbon": 107.0,
            "solubility_sens": 0.03,
            "t_half": 0.2,
        },
    )
    print(decreased_fnpp_results["Biosphere carbon flux"])
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_carbon={"ml_fracmax": 0.7},
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
            "mixed_carbon": 107.0,
            "ml_fracmax": 0.7,
            "solubility_sens": 0.03,
            "t_half": 0.2,
        },
    )
    # TODO check these to make sure the results and checks make sense
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

    # Plots section that you can comment in to look at plots

    """
    fig, axs = plt.subplots(nrows = 2, ncols=3, sharex=True, figsize=(16,8))
    test_data = {
        "default": default_results,
        "beta_f": beta_f_results,
        "mixed": mixed_carbon_results,
        "increasing_fnpp": increased_fnpp_results,
        "decreasing_fnpp": decreased_fnpp_results
        }
    for name, data_dict in test_data.items():
        for i, (var_name, data) in enumerate(data_dict.items()):
            axs[i//3, i%3].plot(np.arange(nystart, nyend+1), data, label = name)
            axs[i//3, i%3].set_title(var_name)
            axs[i//3, i%3].set_xlabel("Year")

    axs[1, 1].legend()
    fig.savefig("test_plot_full_ts.png")
    plt.clf()

    fig2, axs2 = plt.subplots(nrows = 2, ncols=3, sharex=True, figsize=(16,8))
    for name, data_dict in test_data.items():
        if name in ["beta_f", "mixed", "default"]:
            continue
        print(data_dict.keys())
        for i, (var_name, data) in enumerate(data_dict.items()):
            axs2[i//3, i%3].plot(np.arange(nystart, nyend+1), (data-default_results[var_name])/default_results[var_name], label = name)

            #axs[i//2, i%2].plot(np.arange(nystart, nyend+1), (data-default_results[var_name])/default_results[var_name], label = name)
            axs2[i//3, i%3].set_title(var_name)
            axs2[i//3, i%3].set_xlabel("Year")

    axs2[1, 1].legend()
    fig2.savefig("test_plot_rel_error.png")
    """
    # assert np.all(
    #    increased_fnpp_results["Temperature"] < decreased_fnpp_results["Temperature"]
    # )
    # assert np.all(
    #    default_results["Temperature"] < decreased_fnpp_results["Temperature"]
    # )

    # TODO: Slightly more meaningful tests for these or other temperature varrying changes
    for key, value in increased_fnpp_results.items():
        assert not np.allclose(value, default_results[key])

    for key, value in decreased_fnpp_results.items():
        if key == "CO2 concentration":
            continue
        assert not np.allclose(value, default_results[key])
