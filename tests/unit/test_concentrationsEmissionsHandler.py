import os

from ciceroscm import concentrations_emissions_handler


def test_read_components(test_data_dir):
    df_gas = concentrations_emissions_handler.read_components(
        os.path.join(test_data_dir, "gases_v1RCMIP.txt")
    )
    print(df_gas.columns)
    assert df_gas.columns.tolist() == [
        "EM_UNIT",
        "CONC_UNIT",
        "BETA",
        "ALPHA",
        "TAU1",
        "TAU2",
        "TAU3",
        "NAT_EM",
    ]
    assert len(df_gas.index) == 46


def test_rs_and_rb_functions():
    assert concentrations_emissions_handler._rs_function(0) == 1.0
    assert concentrations_emissions_handler._rs_function(3.0) == 0.6884435390654896
    assert concentrations_emissions_handler._rb_function(0) == -3.599999999991197e-06
