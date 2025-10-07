import os

import numpy as np
import pandas as pd

from ciceroscm import CICEROSCM, input_handler
from ciceroscm.carbon_cycle.carbon_cycle_mod import (
    CarbonCycleModel as DefaultCarbonCycleModel,
)
from ciceroscm.carbon_cycle.carbon_cycle_mod_box import (
    CarbonCycleModel as BoxCarbonCycleModel,
)
from ciceroscm.carbon_cycle.carmen_wrapper import (
    CarbonCycleModel as CarmenCarbonCycleModelWrapper
)
from ciceroscm.carbon_cycle.carmen.src.carmen import (
    CarbonCycle as CarmenCarbonCycleModel
)
from ciceroscm.two_layer_ocean import TwoLayerOceanModel
from ciceroscm.upwelling_diffusion_model import UpwellingDiffusionModel


import matplotlib.pyplot as plt

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


def test_ciceroscm_with_carmen_carbon_cycle_model(test_data_dir):

    # co2_emissions = pd.read_csv("/home/eearp/code/coupling-carmen-cscm/ciceroscm/src/ciceroscm/carbon_cycle/carmen/src/carmen/data/scenarios/sce_UKESM1-0-LL_ssp245.txt", sep="\s+")["emission"]

    ih = input_handler.InputHandler({"nystart": 1750, "nyend": 2100, "emstart": 1850})
    emissions_data = ih.read_emissions(
        os.path.join(test_data_dir, "ssp245_em_RCMIP.txt")
    )

    # emissions_data.iloc[100:, 0] = co2_emissions.values
    # emissions_data.iloc[:, 1] = np.zeros(len(emissions_data.index))

    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            # "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_data": emissions_data,
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
            "carbon_cycle_model": "carmen",  # Specify the box model
            "conc_run": False,
        },
    )
    assert isinstance(cscm.ce_handler.carbon_cycle, CarmenCarbonCycleModelWrapper)
    assert isinstance(cscm.ce_handler.carbon_cycle.carmen, CarmenCarbonCycleModel)
    cscm._run({"results_as_dict": True, "carbon_cycle_outputs": True})
    # print(cscm.results["carbon cycle"])
    assert isinstance(cscm.results["carbon cycle"], pd.DataFrame)
    assert "Biosphere carbon flux" in cscm.results["carbon cycle"].columns
    assert "Biosphere carbon pool" in cscm.results["carbon cycle"].columns
    assert "Vegetation carbon pool" in cscm.results["carbon cycle"].columns
    assert "Soil carbon pool" in cscm.results["carbon cycle"].columns
    assert "Ocean carbon pool" in cscm.results["carbon cycle"].columns
    assert "Net primary production" in cscm.results["carbon cycle"].columns
    assert "Litterfall" in cscm.results["carbon cycle"].columns
    assert "Soil respiration" in cscm.results["carbon cycle"].columns
    assert "Ocean carbon flux" in cscm.results["carbon cycle"].columns

    # assert "asd" == cscm.results.keys()

    for var in cscm.results["carbon cycle"].columns:
        plt.plot(np.arange(1750, 2101), cscm.results["carbon cycle"][var], label=var)
        plt.legend()
        # plt.xlim(1850, 1900)
        plt.savefig(f"{var}_NORESM.png")
        plt.clf()
    
    plt.plot(np.arange(1750, 2101), cscm.results["concentrations"]["CO2"])
    # plt.xlim(1850, 1900)
    # plt.ylim(0, 500)
    plt.savefig(f"CO2_conc_NORESM.png")
    plt.clf()

    plt.plot(np.arange(1750, 2101), cscm.results["emissions"]["CO2"])
    # plt.xlim(1850, 1900)
    # plt.ylim(0, 500)
    plt.savefig(f"CO2_emis_NORESM.png")
    plt.clf()

    plt.plot(np.arange(1750, 2101), cscm.results["dT_glob"])
    # plt.xlim(1850, 1900)
    # plt.ylim(0, 500)
    plt.savefig(f"temperature_NORESM.png")
    plt.clf()


    cscm_vanilla = CICEROSCM(
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
            # "carbon_cycle_model": "box",  # Specify the box model
            "conc_run": False,
        },
    )
    cscm_vanilla._run({"results_as_dict": True, "carbon_cycle_outputs": True})

    for var in cscm_vanilla.results["carbon cycle"].columns:
        plt.plot(np.arange(1750, 2101), cscm_vanilla.results["carbon cycle"][var], label=var)
        plt.legend()
        # plt.xlim(1850, 1900)
        plt.savefig(f"{var}_vanilla_NORESM.png")
        plt.clf()
    
    plt.plot(np.arange(1750, 2101), cscm_vanilla.results["concentrations"]["CO2"])
    # plt.xlim(1850, 1900)
    # plt.ylim(0, 500)
    plt.savefig(f"CO2_conc_vanilla_NORESM.png")
    plt.clf()

    plt.plot(np.arange(1750, 2101), cscm_vanilla.results["emissions"]["CO2"])
    # plt.xlim(1850, 1900)
    # plt.ylim(0, 500)
    plt.savefig(f"CO2_emis_vanilla_NORESM.png")
    plt.clf()

    plt.plot(np.arange(1750, 2101), cscm_vanilla.results["dT_glob"])
    # plt.xlim(1850, 1900)
    # plt.ylim(0, 500)
    plt.savefig(f"temperature_vanilla_NORESM.png")
    plt.clf()

    # print(cscm.results["carbon cycle"]["Biosphere carbon pool"])
    # assert [0] == cscm.results["carbon cycle"]["Net primary production"][-1]
    # back_calc_npp = cscm.ce_handler.carbon_cycle.calculate_npp(
    #     cscm.results["dT_glob"],
    #     co2_conc_series=cscm.results["concentrations"]["CO2"].values,
    # )
    # back_calc_ocean = cscm.ce_handler.carbon_cycle.calculate_ocean_uptake(
    #     cscm.results["dT_glob"],
    #     co2_conc_series=cscm.results["concentrations"]["CO2"].values,
    # )
    # assert np.allclose(
    #     back_calc_npp, cscm.results["carbon cycle"]["Net primary production"].values
    # )
    # assert np.allclose(
    #     back_calc_ocean, cscm.results["carbon cycle"]["Ocean carbon flux"].values
    # )


def test_ciceroscm_with_default_thermal_model(test_data_dir):
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
    assert isinstance(cscm.thermal(), UpwellingDiffusionModel)


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
    assert isinstance(cscm.thermal(), TwoLayerOceanModel)
    cscm._run({"results_as_dict": True})
    assert len(cscm.results["RIB_glob"]) == 351
