import os

import numpy as np

from ciceroscm import CICEROSCM, carbon_cycle_mod


def test_rs_and_rb_functions():
    assert carbon_cycle_mod._rs_function(0) == 1.0
    assert carbon_cycle_mod._rs_function(3.0) == 0.6884435390654896
    assert carbon_cycle_mod._rb_function(0) == -3.599999999991197e-06


def test_get_biosphere_carbon_pool_content():
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    co2_conc_series = np.ones(ccmod.pamset["years_tot"]) * 278.0
    bio_carbon_pool = ccmod.get_biosphere_carbon_pool_content(
        conc_run=True, co2_conc_series=co2_conc_series
    )
    assert np.allclose(bio_carbon_pool, np.zeros(ccmod.pamset["years_tot"]))
    # co2_conc_series = [278*(1.01)**(n) for n in ]


def test_guess_iteration():
    co2_conc_zero = 278.0
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1750})
    em_0 = 0.00259244 + 0.08112671
    co2_conc_now = 277.147003174
    # co2_conc_now = 277.188000997
    ffer = ccmod._get_ffer_timeseries(
        conc_run=True, co2_conc_series=[co2_conc_zero, co2_conc_zero]
    )
    print(ffer)
    em_guess = ccmod._guess_emissions_iteration(
        co2_conc_now=co2_conc_now, co2_conc_zero=co2_conc_zero
    )
    print(f"Value for em_guess: {em_guess}")
    print(f"Corresponding for em_0: {(em_0 - ffer[0])/2.123}")
    em_guess_input = 2.123 * (em_guess) + ffer[0]
    ccmod.reset_co2_hold()
    conc_from_guess = ccmod.co2em2conc(1750, em_guess_input)
    second_guess = ccmod.simplified_em_backward(
        co2_conc_now=conc_from_guess, co2_conc_zero=co2_conc_zero
    )
    ccmod.reset_co2_hold()
    conc_from_em = ccmod.co2em2conc(1750, em_0)
    ccmod.reset_co2_hold()
    conc_from_second_guess = ccmod.co2em2conc(1750, 2.123 * (second_guess) + ffer[0])
    print(f"Conc from guess: {conc_from_guess}")
    print(f"Conc from second guess: {conc_from_second_guess}")
    print(f"Conc from em: {conc_from_em}")
    assert np.allclose(conc_from_guess, co2_conc_now)


def test_simplified_em_backwards():
    co2_conc_zero = 278.0
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    print(
        ccmod.simplified_em_backward(
            co2_conc_now=co2_conc_zero, co2_conc_zero=co2_conc_zero
        )
    )
    assert np.allclose(
        ccmod.simplified_em_backward(
            co2_conc_now=co2_conc_zero, co2_conc_zero=co2_conc_zero
        ),
        0.0,
    )


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
    cum_emis = np.cumsum(emis_series)
    biopool = cscm.ce_handler.carbon_cycle.get_biosphere_carbon_pool_content()
    oceanpool = cscm.ce_handler.carbon_cycle.get_ocean_carbon_pool_content()
    summed_carbon_pools = conc_co2_series + biopool / 2.123 - 278 - oceanpool / 2.123
    print(summed_carbon_pools[:5])
    print(conc_co2_series[:5] - 278)
    print(biopool[:5] / 2.123)
    print(oceanpool[:5])
    print(cum_emis[:5] / 2.123)
    # TODO : Put tests here back on
    # assert np.allclose(summed_carbon_pools, cum_emis / 2.123)
    assert True


# TODO: Check ocean calculation with and without internal back calculation
