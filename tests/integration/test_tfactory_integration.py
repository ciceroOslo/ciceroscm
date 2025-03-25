import pytest
import os
from ciceroscm import CICEROSCM
from ciceroscm.upwelling_diffusion_model import UpwellingDiffusionModel
from ciceroscm.two_layer_ocean import TwoLayerOceanModel


def test_ciceroscm_with_default_thermal_model(tmpdir, test_data_dir):
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
            "thermal_model": "default",  # Specify the default model
        },  # lambda it, idtm: 0.5 + 0.5 * np.exp(-it / idtm / 100.0),
    )
    assert isinstance(cscm.thermal({}), UpwellingDiffusionModel)


def test_ciceroscm_with_twolayer_thermal_model(tmpdir, test_data_dir):
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
            "thermal_model": "twolayer",  # Specify the box model
        },
    )
    assert isinstance(cscm.thermal({}), TwoLayerOceanModel)
