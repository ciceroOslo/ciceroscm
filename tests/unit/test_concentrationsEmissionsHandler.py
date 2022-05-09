from ciceroscm import concentrations_emissions_handler


def test_rs_and_rb_functions():
    assert concentrations_emissions_handler._rs_function(0) == 1.0
    assert concentrations_emissions_handler._rs_function(3.0) == 0.6884435390654896
    assert concentrations_emissions_handler._rb_function(0) == -3.599999999991197e-06
