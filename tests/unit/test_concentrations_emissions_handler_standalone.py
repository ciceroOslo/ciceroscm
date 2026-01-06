import os

from ciceroscm import concentrations_emissions_handler, input_handler


def test_ce_handler_init(test_data_dir):
    ce_handler = concentrations_emissions_handler.ConcentrationsEmissionsHandler(
        input_handler.InputHandler(
            {
                "gaspam_file": os.path.join(
                    test_data_dir, "gases_vupdate_2022_AR6.txt"
                ),
                "concentrations_file": os.path.join(
                    test_data_dir, "ssp245_conc_RCMIP.txt"
                ),
                "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
                "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
                "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
                "nystart": 1750,
                "nyend": 2100,
            }
        ),
        pamset={"lifetime_mode": "CONSTANT_12"},
    )
    assert ce_handler.pamset["carbon_cycle_model"] == "default"
    ce_handler.reset_with_new_pams({"qo3": 0.4})
    assert ce_handler.pamset["qo3"] == 0.4
    assert ce_handler.pamset["qnh3"] == 0.0
