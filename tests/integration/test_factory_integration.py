import os

import numpy as np
import pandas as pd

from ciceroscm import CICEROSCM
from ciceroscm.carbon_cycle.carbon_cycle_mod import (
    CarbonCycleModel as DefaultCarbonCycleModel,
)
from ciceroscm.carbon_cycle.carbon_cycle_mod_box import (
    CarbonCycleModel as BoxCarbonCycleModel,
)
from ciceroscm.thermal_model.two_layer_ocean import TwoLayerOceanModel
from ciceroscm.thermal_model.upwelling_diffusion_model import UpwellingDiffusionModel


def test_ciceroscm_with_default_carbon_cycle_model(test_data_dir):
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


def test_ciceroscm_with_box_carbon_cycle_model(test_data_dir):
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
            "conc_run": False,
        },
    )
    assert isinstance(cscm.ce_handler.carbon_cycle, BoxCarbonCycleModel)
    cscm._run({"results_as_dict": True, "carbon_cycle_outputs": True})
    # print(cscm.results["carbon cycle"])
    assert isinstance(cscm.results["carbon cycle"], pd.DataFrame)
    assert "Ocean carbon flux" in cscm.results["carbon cycle"].columns
    back_calc_npp = cscm.ce_handler.carbon_cycle.calculate_npp(
        cscm.results["dT_glob"],
        co2_conc_series=cscm.results["concentrations"]["CO2"].values,
    )
    back_calc_ocean = cscm.ce_handler.carbon_cycle.calculate_ocean_uptake(
        cscm.results["dT_glob"],
        co2_conc_series=cscm.results["concentrations"]["CO2"].values,
    )
    assert np.allclose(
        back_calc_npp, cscm.results["carbon cycle"]["Net primary production"].values
    )
    assert np.allclose(
        back_calc_ocean, cscm.results["carbon cycle"]["Ocean carbon flux"].values
    )


def test_ciceroscm_with_default_thermal_model(test_data_dir):
    """
    Test that CICEROSCM uses the default thermal model when configured.
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
    # Verify the thermal model configuration is stored correctly
    assert cscm.cfg["thermal_model"] == "default"

    # Verify the thermal model is created correctly during run
    from ciceroscm.component_factory_functions import create_thermal_model

    thermal_model_class = create_thermal_model(cscm.cfg["thermal_model"])
    assert thermal_model_class == UpwellingDiffusionModel


def test_ciceroscm_with_twolayer_thermal_model(test_data_dir):
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
    # Verify the thermal model configuration is stored correctly
    assert cscm.cfg["thermal_model"] == "twolayer"

    # Verify the thermal model is created correctly during factory call
    from ciceroscm.component_factory_functions import create_thermal_model

    thermal_model_class = create_thermal_model(cscm.cfg["thermal_model"])
    assert thermal_model_class == TwoLayerOceanModel

    # Run the model to ensure it works with the two-layer thermal model
    cscm._run({"results_as_dict": True})
    assert len(cscm.results["RIB_glob"]) == 351
