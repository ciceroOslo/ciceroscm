import os

from ciceroscm import concentrations_emissions_handler, perturbations


def test_check_numeric_pamset(test_data_dir):
    emissions_file = os.path.join(test_data_dir, "ssp245_em_RCMIP.txt")
    pert_file = os.path.join(test_data_dir, "pertem_test.txt")
    em_df = concentrations_emissions_handler.read_inputfile(emissions_file)
    em_df_change = perturbations.perturb_emissions(pert_file, em_df)
    print(em_df["CO2"][2000])
    print(em_df_change["CO2"][2000])
    em_df_old = concentrations_emissions_handler.read_inputfile(emissions_file)
    diff = em_df_old.sub(em_df_change)
    print(diff["CO2"][2000])
    assert diff["CO2"][1999] == 0.0
    assert diff["CO2"][2000] == 4.0
    assert diff["CO2"][2001] == 4.0
    assert diff["CO2"][2002] == 4.0
    assert diff["CH4"][1864] == 0.0
    assert diff["CH4"][2000] == 4.0
    assert diff["CH4"][2001] == 4.0
    assert diff["CH4"][2002] == 5.0
    assert diff["N2O"][2078] == 0.0
    assert diff["N2O"][2000] == -8.0
    assert diff["N2O"][2001] == 4.0
    assert diff["N2O"][2002] == 4.0
