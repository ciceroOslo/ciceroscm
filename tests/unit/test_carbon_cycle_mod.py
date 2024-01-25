from ciceroscm import carbon_cycle_mod


def test_rs_and_rb_functions():
    assert carbon_cycle_mod._rs_function(0) == 1.0
    assert carbon_cycle_mod._rs_function(3.0) == 0.6884435390654896
    assert carbon_cycle_mod._rb_function(0) == -3.599999999991197e-06
