import os

import numpy as np

from ciceroscm import CICEROSCM
from ciceroscm.carbon_cycle import carbon_cycle_mod


def test_linear_fnpp_from_temp():
    assert carbon_cycle_mod.linear_fnpp_from_temp() == 60.0
    assert carbon_cycle_mod.linear_fnpp_from_temp(fnpp_temp_coeff=1) == 60.0
    assert carbon_cycle_mod.linear_fnpp_from_temp(dtemp=1) == 60.0
    assert carbon_cycle_mod.linear_fnpp_from_temp(fnpp_temp_coeff=2, dtemp=3) == 66.0


def test_default_pamset_values(test_data_dir):
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    assert ccmod.pamset["beta_f"] == 0.287
    assert ccmod.pamset["mixed_carbon"] == 75.0
    assert ccmod.pamset["fnpp_temp_coeff"] == 0
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
    assert ccmod_inside.pamset["fnpp_temp_coeff"] == 0


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
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1750})
    co2_conc_now = 277.147003174
    em_guess = ccmod._guess_emissions_iteration(
        co2_conc_now=co2_conc_now, co2_conc_zero=co2_conc_zero
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
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2100, "nystart": 1750})
    em_back_calculated = ccmod.back_calculate_emissions(conc_co2_series)
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
    conc_co2_series_no_feedback = cscm.results["concentrations"]["CO2"].values
    print(cscm.ce_handler.carbon_cycle.pamset)
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_carbon={"fnpp_temp_coeff": -10},
    )
    conc_co2_series = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    temp_timseries = cscm.results["dT_glob"]
    print(cscm.ce_handler.carbon_cycle.pamset)

    ccmod = carbon_cycle_mod.CarbonCycleModel(
        {"nyend": 2100, "nystart": 1750}, pamset_carbon={"fnpp_temp_coeff": -10}
    )
    em_back_calculated = ccmod.back_calculate_emissions(
        conc_co2_series, dtemp_series=temp_timseries
    )
    assert not np.allclose(conc_co2_series, conc_co2_series_no_feedback)
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
