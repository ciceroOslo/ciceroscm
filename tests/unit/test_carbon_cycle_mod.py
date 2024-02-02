import numpy as np

from ciceroscm import carbon_cycle_mod


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
