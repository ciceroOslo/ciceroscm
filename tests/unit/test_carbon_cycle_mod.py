import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from ciceroscm import CICEROSCM
from ciceroscm.carbon_cycle import carbon_cycle_mod
from ciceroscm import input_handler


def test_linear_fnpp_from_temp():
    assert carbon_cycle_mod.linear_fnpp_from_temp() == 60.0
    assert carbon_cycle_mod.linear_fnpp_from_temp(fnpp_temp_coeff=1) == 60.0
    assert carbon_cycle_mod.linear_fnpp_from_temp(dtemp=1) == 60.0
    assert carbon_cycle_mod.linear_fnpp_from_temp(fnpp_temp_coeff=2, dtemp=3) == 66.0


def test_default_pamset_values(test_data_dir):
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    assert ccmod.pamset["beta_f"] == 0.287
    assert ccmod.pamset["mixed_carbon"] == 75.0
    assert ccmod.pamset["fnpp_temp_coeff"] == 0
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
    ccmod_inside = cscm.ce_handler.carbon_cycle
    assert ccmod_inside.pamset["beta_f"] == 0.287
    assert ccmod_inside.pamset["mixed_carbon"] == 75.0
    assert ccmod_inside.pamset["fnpp_temp_coeff"] == 0


def test_get_biosphere_carbon_flux():
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1850})
    co2_conc_series = (
        np.ones(ccmod.pamset["years_tot"]) * carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC
    )
    bio_carbon_flux = ccmod.get_biosphere_carbon_flux(
        conc_run=True, co2_conc_series=co2_conc_series
    )
    assert np.allclose(bio_carbon_flux, np.zeros(ccmod.pamset["years_tot"]))
    # co2_conc_series = [278*(1.01)**(n) for n in ]


def test_guess_iteration():
    co2_conc_zero = carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC
    ccmod = carbon_cycle_mod.CarbonCycleModel({"nyend": 2015, "nystart": 1750})
    co2_conc_now = 277.147003174
    em_guess = ccmod._guess_emissions_iteration(
        co2_conc_now=co2_conc_now, co2_conc_zero=co2_conc_zero
    )
    print(f"Value for em_guess: {em_guess}")
    ccmod.reset_co2_hold()
    conc_from_guess = ccmod.co2em2conc(1750, em_guess)
    print(f"Conc from guess: {conc_from_guess}")
    print(f"Conc from em: {co2_conc_now}")
    assert np.allclose(conc_from_guess, co2_conc_now)


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


def test_back_calculate_emissions_with_temperature_feedback(test_data_dir):
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

    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
    )
    conc_co2_series_no_feedback = cscm.results["concentrations"]["CO2"].values
    print(cscm.ce_handler.carbon_cycle.pamset)
    cscm._run(
        {"results_as_dict": True, "carbon_cycle_outputs": True},
        pamset_carbon={"fnpp_temp_coeff": -10},
    )
    conc_co2_series = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    temp_timseries = cscm.results["dT_glob"]
    print(cscm.ce_handler.carbon_cycle.pamset)

    ccmod = carbon_cycle_mod.CarbonCycleModel(
        {"nyend": 2100, "nystart": 1750}, pamset_carbon={"fnpp_temp_coeff": -10}
    )
    em_back_calculated = ccmod.back_calculate_emissions(
        conc_co2_series, dtemp_series=temp_timseries
    )
    assert not np.allclose(conc_co2_series, conc_co2_series_no_feedback)
    assert np.allclose(em_back_calculated, emis_series, rtol=1.0e-2)
    # TODO: Test carbon cycle outputs with feedbacks


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
    cscm._run({"results_as_dict": True, "carbon_cycle_outputs": True})
    conc_co2_series = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    cum_emis = np.cumsum(emis_series)
    bioflux = cscm.ce_handler.carbon_cycle.get_biosphere_carbon_flux()
    oceanflux = cscm.ce_handler.carbon_cycle.get_ocean_carbon_flux()
    assert np.allclose(bioflux, cscm.results["carbon cycle"]["Biosphere carbon flux"])
    assert np.allclose(oceanflux, cscm.results["carbon cycle"]["Ocean carbon flux"])
    summed_carbon_pools = (
        (conc_co2_series
        - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)
        + np.cumsum(bioflux) / carbon_cycle_mod.PPM_CO2_TO_PG_C
        + np.cumsum(oceanflux) / carbon_cycle_mod.PPM_CO2_TO_PG_C
    )
    atmospheric_flux = reverse_cumsum((conc_co2_series - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)) * carbon_cycle_mod.PPM_CO2_TO_PG_C
    summed_fluxes = atmospheric_flux + bioflux + oceanflux
    print(summed_carbon_pools[:5])
    print(conc_co2_series[:5] - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)
    print(bioflux[:5] / carbon_cycle_mod.PPM_CO2_TO_PG_C)
    print(oceanflux[:5])
    print(cum_emis[:5] / carbon_cycle_mod.PPM_CO2_TO_PG_C)
    cumulative_mismatch1 = summed_carbon_pools * carbon_cycle_mod.PPM_CO2_TO_PG_C - cum_emis
    cumulative_mismatch2 = -np.cumsum(emis_series - summed_fluxes)
    
    # TODO : Put tests here back on
    #yearly_total_flux = reverse_cumsum()

    fig1, axs = plt.subplots(nrows = 3, ncols = 5)
    years = cscm.results["concentrations"].index
    axs[1,0].plot(years, emis_series, label = "Emissions")
    axs[0,0].plot(years, cum_emis, label = "Cumulative emissions")
    axs[0,1].plot(years, summed_carbon_pools * carbon_cycle_mod.PPM_CO2_TO_PG_C, label = "Cumulative pools")
    axs[0,2].plot(years, cumulative_mismatch1, label = "Cumulative mismatch")
    axs[0,3].plot(years, cumulative_mismatch1 / cum_emis, label = "Cumulative mismatch rel")
    axs[0,4].plot(years[50:], (cumulative_mismatch1 / cum_emis)[50:], label = "Cumulative mismatch rel")
    axs[1,1].plot(years, summed_fluxes, label="Carbon cycle fluxes")
    axs[1,2].plot(years, (summed_fluxes - emis_series), label= "Flux mismatch")
    axs[1,3].plot(years, (summed_fluxes - emis_series) / emis_series, label= "Flux mismatch rel")
    axs[1,4].plot(years[50:], (summed_fluxes - emis_series)[50:] / emis_series[50:], label= "Flux mismatch rel")
    axs[2,0].plot(years, bioflux, label="Biosphere carbon flux")
    axs[2,1].plot(years, oceanflux, label="Ocean carbon flux")
    axs[2,2].plot(years, atmospheric_flux, label="Atmospheric carbon flux")
    axs[2,3].plot(years, cumulative_mismatch1 / emis_series, label = "Cumulative mismatch rel to emissions")
    axs[2,4].plot(years[50:], (cumulative_mismatch1 / emis_series)[50:], label = "Cumulative mismatch rel to emissions")
    for ax in axs.flatten():
        ax.legend()
        ax.set_xlabel("Years")
    fig1.savefig("carbon_budget_plot.png")
    fig2, ax = plt.subplots()
    ax.plot(years, summed_fluxes - emis_series, label= "Flux mismatch")
    ax.plot(years, bioflux, label="Biosphere carbon flux")
    ax.plot(years, oceanflux, label="Ocean carbon flux")
    ax.plot(years, atmospheric_flux, label="Atmospheric carbon flux")
    ax.plot(years, cumulative_mismatch1, label="Cumulative mismatch")
    ax.legend()
    ax.set_xlabel("Years")
    fig2.savefig("all_fluxes_in_one.png")
    print(bioflux[0])
    print(oceanflux[0])
    print(atmospheric_flux[0])
    print((summed_fluxes - emis_series)[0])
    assert np.allclose(cumulative_mismatch2, cumulative_mismatch1)
    assert np.allclose(summed_fluxes, emis_series, rtol = 1e-2)
    assert np.allclose(summed_carbon_pools, (cum_emis / carbon_cycle_mod.PPM_CO2_TO_PG_C))
    assert True

def test_carbon_pools_flat(test_data_dir):
    ih = input_handler.InputHandler(
        {
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
        }
    )
    em_data = (ih.get_data("emissions"))
    em_data_flat_co2 = pd.DataFrame(
        data = np.zeros_like(em_data.values),
        columns=em_data.columns,
        index=em_data.index
    )
    em_data_flat_co2.loc[:, "CO2_FF"] = 5*np.ones(len(em_data.index))
    print(em_data_flat_co2.head())

    cscm = CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_vupdate_2022_AR6.txt"),
            "nyend": 2100,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_data": em_data_flat_co2,
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        },
    )
    cscm._run({"results_as_dict": True, "carbon_cycle_outputs": True})
    conc_co2_series = cscm.results["concentrations"]["CO2"].values
    emis_series = cscm.results["emissions"]["CO2"].values
    cum_emis = np.cumsum(emis_series)
    bioflux = cscm.ce_handler.carbon_cycle.get_biosphere_carbon_flux()
    oceanflux = cscm.ce_handler.carbon_cycle.get_ocean_carbon_flux()
    assert np.allclose(bioflux, cscm.results["carbon cycle"]["Biosphere carbon flux"])
    assert np.allclose(oceanflux, cscm.results["carbon cycle"]["Ocean carbon flux"])
    summed_carbon_pools = (
        (conc_co2_series
        - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)
        + np.cumsum(bioflux) / carbon_cycle_mod.PPM_CO2_TO_PG_C
        + np.cumsum(oceanflux) / carbon_cycle_mod.PPM_CO2_TO_PG_C
    )
    atmospheric_flux = reverse_cumsum((conc_co2_series - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)) * carbon_cycle_mod.PPM_CO2_TO_PG_C
    summed_fluxes = atmospheric_flux + bioflux + oceanflux
    print(summed_carbon_pools[:5])
    print(conc_co2_series[:5] - carbon_cycle_mod.PREINDUSTRIAL_CO2_CONC)
    print(bioflux[:5] / carbon_cycle_mod.PPM_CO2_TO_PG_C)
    print(oceanflux[:5])
    print(cum_emis[:5] / carbon_cycle_mod.PPM_CO2_TO_PG_C)
    cumulative_mismatch1 = summed_carbon_pools * carbon_cycle_mod.PPM_CO2_TO_PG_C - cum_emis
    cumulative_mismatch2 = -np.cumsum(emis_series - summed_fluxes)
    
    # TODO : Put tests here back on
    #yearly_total_flux = reverse_cumsum()

    fig1, axs = plt.subplots(nrows = 3, ncols = 3)
    years = cscm.results["concentrations"].index
    axs[1,0].plot(years, emis_series, label = "Emissions")
    axs[0,0].plot(years, cum_emis, label = "Cumulative emissions")
    axs[0,1].plot(years, summed_carbon_pools * carbon_cycle_mod.PPM_CO2_TO_PG_C, label = "Cumulative pools")
    axs[0,2].plot(years, cumulative_mismatch1, label = "Cumulative mismatch")
    #axs[0,3].plot(years, cumulative_mismatch1 / cum_emis, label = "Cumulative mismatch rel")
    #axs[0,4].plot(years[100:], (cumulative_mismatch1 / cum_emis)[100:], label = "Cumulative mismatch rel")
    axs[1,1].plot(years, summed_fluxes, label="Carbon cycle fluxes")
    axs[1,2].plot(years[:], (summed_fluxes - emis_series)[:], label= "Flux mismatch")
    #axs[1,3].plot(years, (summed_fluxes - emis_series) / emis_series, label= "Flux mismatch rel")
    #axs[1,4].plot(years[100:], (summed_fluxes - emis_series)[100:] / emis_series[50:], label= "Flux mismatch rel")
    axs[2,0].plot(years, bioflux, label="Biosphere carbon flux")
    axs[2,1].plot(years, oceanflux, label="Ocean carbon flux")
    axs[2,2].plot(years, atmospheric_flux, label="Atmospheric carbon flux")
    #axs[2,3].plot(years, cumulative_mismatch1 / emis_series, label = "Cumulative mismatch rel to emissions")
    #axs[2,4].plot(years[50:], (cumulative_mismatch1 / emis_series)[50:], label = "Cumulative mismatch rel to emissions")
    for ax in axs.flatten():
        ax.legend()
        ax.set_xlabel("Years")
    fig1.savefig("carbon_budget_plot_flat.png")
    fig2, ax = plt.subplots()
    ax.plot(years, summed_fluxes - emis_series, label= "Flux mismatch")
    ax.plot(years, bioflux, label="Biosphere carbon flux")
    ax.plot(years, oceanflux, label="Ocean carbon flux")
    ax.plot(years, atmospheric_flux, label="Atmospheric carbon flux")
    ax.plot(years, cumulative_mismatch1, label="Cumulative mismatch")
    ax.legend()
    ax.set_xlabel("Years")
    fig2.savefig("all_fluxes_in_one_flat.png")
    print(bioflux[0])
    print(oceanflux[0])
    print(atmospheric_flux[0])
    print((summed_fluxes - emis_series)[0])
    assert np.allclose(cumulative_mismatch2, cumulative_mismatch1)
    assert np.allclose(summed_fluxes, emis_series, rtol = 1e-2)
    assert np.allclose(summed_carbon_pools, (cum_emis / carbon_cycle_mod.PPM_CO2_TO_PG_C))
    assert True


def reverse_cumsum(cumulated):
    print(cumulated[:-1].copy())
    cumsum_shifted = np.insert(cumulated[:-1].copy(), 0, 0)
    print(cumulated[0])
    decumulated = cumulated - cumsum_shifted
    print(decumulated[0])
    return decumulated

# TODO: Check ocean calculation with and without internal back calculation
