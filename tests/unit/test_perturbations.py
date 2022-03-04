import os

import numpy as np

from ciceroscm import concentrations_emissions_handler, perturbations


def test_emission_perturbation(test_data_dir):
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


def test_make_forcing_perturbations_dataset(test_data_dir):
    pert_file = os.path.join(test_data_dir, "pertforc_test.txt")
    forc_pert = perturbations.ForcingPerturbation(pert_file, 1999)
    assert forc_pert.check_if_year_in_pert(2000)
    assert not forc_pert.check_if_year_in_pert(1983)
    assert forc_pert.check_if_compound_in_pert("OTHER")
    assert not forc_pert.check_if_compound_in_pert("CO2")
    test = np.array([38.0, 42.0, 48.0, 48.0])
    forc = {"CO2": test, "OTHER": [0, 0, 0, 0]}

    for year in range(1999, 2003):
        totforc = forc["CO2"][year - 1999]
        totforc, forc_nh, forc_sh, forc = forc_pert.add_forcing_pert(
            totforc, totforc, totforc, forc, year
        )
        if year == 1999:
            expected = test[0]
        else:
            expected = test[year - 1999] - 4.0
        assert totforc == expected
    assert np.array_equal(forc["OTHER"], [0, -4.0, -4.0, -4.0])
