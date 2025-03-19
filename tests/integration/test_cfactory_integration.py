import os

from ciceroscm import CICEROSCM
from ciceroscm.carbon_cycle_mod import CarbonCycleModel as DefaultCarbonCycleModel
from ciceroscm.carbon_cycle_mod_box import CarbonCycleModel as BoxCarbonCycleModel


def test_ciceroscm_with_default_carbon_cycle_model(tmpdir, test_data_dir):
    """
    Test that CICEROSCM uses the default CarbonCycleModel when configured.
    """
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            # "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "carbon_cycle_model": "default",  # Specify the default model
        },  # lambda it, idtm: 0.5 + 0.5 * np.exp(-it / idtm / 100.0),
    )
    assert isinstance(cscm.ce_handler.carbon_cycle, DefaultCarbonCycleModel)


def test_ciceroscm_with_box_carbon_cycle_model(tmpdir, test_data_dir):
    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            # "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "carbon_cycle_model": "box",  # Specify the box model
        },
    )
    assert isinstance(cscm.ce_handler.carbon_cycle, BoxCarbonCycleModel)
