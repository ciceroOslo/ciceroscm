import os

import numpy as np

from ciceroscm import CICEROSCM
from ciceroscm.carbon_cycle import carbon_cycle_mod


def test_default_pamset_values(test_data_dir):
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    assert ccmod.pamset["beta_f"] == 0.287
    assert ccmod.pamset["mixed_carbon"] == 75.0
    assert ccmod.pamset["ml_w_sigmoid"] == 3.0
    assert ccmod.pamset["ml_fracmax"] == 0.5
    assert ccmod.pamset["ml_t_half"] == 0.5
    assert ccmod.pamset["t_half"] == 0.5
    assert ccmod.pamset["w_sigmoid"] == 7
    assert ccmod.pamset["t_threshold"] == 4
    assert ccmod.pamset["w_threshold"] == 7
    assert ccmod.pamset["solubility_sens"] == 0.02
    assert ccmod.pamset["solubility_limit"] == 0.5

    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    ccmod_inside = cscm.ce_handler.carbon_cycle
    assert ccmod_inside.pamset["beta_f"] == 0.287
    assert ccmod_inside.pamset["mixed_carbon"] == 75.0
    assert ccmod_inside.pamset["ml_w_sigmoid"] == 3.0
    assert ccmod_inside.pamset["ml_fracmax"] == 0.5
    assert ccmod_inside.pamset["ml_t_half"] == 0.5
    assert ccmod_inside.pamset["t_half"] == 0.5
    assert ccmod_inside.pamset["w_sigmoid"] == 7
    assert ccmod_inside.pamset["t_threshold"] == 4
    assert ccmod_inside.pamset["w_threshold"] == 7
    assert ccmod_inside.pamset["solubility_sens"] == 0.02
    assert ccmod_inside.pamset["solubility_limit"] == 0.5


def test_get_biosphere_carbon_flux():
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    co2_conc_series = (
        np.ones(ccmod.pamset["years_tot"]) * carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC
    )
    bio_carbon_flux = ccmod.get_biosphere_carbon_flux(
        conc_run=True, co2_conc_series=co2_conc_series
    )
    assert np.allclose(bio_carbon_flux, np.zeros(ccmod.pamset["years_tot"]))
    # co2_conc_series = [278*(1.01)**(n) for n in ]


def test_guess_iteration():
    co2_conc_zero = carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC
    print(co2_conc_zero)
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1750})
    co2_conc_now = 277.147003174
    initial_guess = ccmod.get_initial_max_min_guess(co2_conc_now, co2_conc_zero)
    print(initial_guess)
    print()
    em_guess = ccmod._guess_emissions_iteration(
        co2_conc_now=co2_conc_now,
        initial_max_min_guess=initial_guess,
    )
    print(f"Value for em_guess: {em_guess}")
    ccmod.reset_co2_hold()
    conc_from_guess = ccmod.co2em2conc(1750, em_guess)
    print(f"Conc from guess: {conc_from_guess}")
    print(f"Conc from em: {co2_conc_now}")
    assert np.allclose(conc_from_guess, co2_conc_now)


def test_back_calculate_emissions(test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    cscm._run({"results_as_dict": True})
    conc_co2_series = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    dtemp_series = cscm.results["dT_glob"]
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2100, "nystart": 1750})
    em_back_calculated = ccmod.back_calculate_emissions(
        conc_co2_series, feedback_dict_series={"dtemp": dtemp_series}
    )
    assert np.allclose(em_back_calculated, emis_series, rtol=1.0e-2)


def test_back_calculate_emissions_with_temperature_feedback(test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )

    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
    )
    conc_co2_series_default = cscm.results["concentrations"]["CO2"].values
    print(cscm.ce_handler.carbon_cycle.pamset)
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_carbon={"t_threshold": 2, "w_threshold": 2},
    )
    conc_co2_series_all_die = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    temp_timseries = cscm.results["dT_glob"]
    print(cscm.ce_handler.carbon_cycle.pamset)

    ccmod = carbon_cycle_mod.CarbonCycleModel(
        {"nyend": 2100, "nystart": 1750},
        pamset_carbon={"t_threshold": 2, "w_threshold": 2},
    )
    em_back_calculated = ccmod.back_calculate_emissions(
        conc_co2_series_all_die, feedback_dict_series={"dtemp": temp_timseries}
    )
    assert not np.allclose(conc_co2_series_all_die, conc_co2_series_default)
    assert np.allclose(em_back_calculated, emis_series, rtol=1.0e-2)
    # TODO: Test carbon cycle outputs with feedbacks


def test_carbon_pools(test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    cscm._run({"results_as_dict": True})
    conc_co2_series = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    bioflux = cscm.ce_handler.carbon_cycle.get_biosphere_carbon_flux()
    oceanflux = cscm.ce_handler.carbon_cycle.get_ocean_carbon_flux()
    atmospheric_flux = (
        reverse_cumsum((conc_co2_series - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC))
        * carbon_cycle_mod.PPM_CO2_TO_PG_C
    )
    summed_fluxes = atmospheric_flux + bioflux + oceanflux
    # print(summed_carbon_pools[:5])
    # print(conc_co2_series[:5] - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)
    # print(bioflux[:5] / carbon_cycle_mod.PPM_CO2_TO_PG_C)
    # print(oceanflux[:5])
    # print(cum_emis[:5] / carbon_cycle_mod.PPM_CO2_TO_PG_C)
    # This is as close we are likely going to get here...
    assert np.allclose(summed_fluxes, emis_series, rtol=1e-2)


def reverse_cumsum(cumulated):
    print(cumulated[:-1].copy())
    cumsum_shifted = np.insert(cumulated[:-1].copy(), 0, 0)
    print(cumulated[0])
    decumulated = cumulated - cumsum_shifted
    print(decumulated[0])
    return decumulated


# TODO: Check ocean calculation with and without internal back calculation
