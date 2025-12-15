"""
Module that reads in CICERO-SCM results
and returns data to append to SCMRun
"""

import numpy as np
import pandas as pd

openscm_to_cscm_dict = {
    "Surface Air Temperature Change": "dT_glob_air",
    # GMST
    "Surface Air Ocean Blended Temperature Change": "dT_glob",
    # ERFs
    "Effective Radiative Forcing": "Total_forcing+sunvolc",
    "Effective Radiative Forcing|Anthropogenic": "Total_forcing",
    "Effective Radiative Forcing|Aerosols": "Aerosols",
    "Effective Radiative Forcing|Aerosols|Direct Effect": "Aerosols|Direct Effect",
    "Effective Radiative Forcing|Aerosols|Direct Effect|BC": "BC",
    "Effective Radiative Forcing|Aerosols|Direct Effect|OC": "OC",
    "Effective Radiative Forcing|Aerosols|Direct Effect|SOx": "SO2",
    "Effective Radiative Forcing|Aerosols|Indirect Effect": "SO4_IND",
    "Effective Radiative Forcing|Greenhouse Gases": "GHG",
    "Effective Radiative Forcing|F-Gases": "Fgas",
    "Effective Radiative Forcing|HFC125": "HFC125",
    "Effective Radiative Forcing|HFC134a": "HFC134a",
    "Effective Radiative Forcing|HFC143a": "HFC143a",
    "Effective Radiative Forcing|HFC227ea": "HFC227ea",
    "Effective Radiative Forcing|HFC23": "HFC23",
    "Effective Radiative Forcing|HFC245fa": "HFC245fa",
    "Effective Radiative Forcing|HFC32": "HFC32",
    "Effective Radiative Forcing|HFC4310mee": "HFC4310mee",
    "Effective Radiative Forcing|CF4": "CF4",
    "Effective Radiative Forcing|C6F14": "C6F14",
    "Effective Radiative Forcing|C2F6": "C2F6",
    "Effective Radiative Forcing|SF6": "SF6",
    "Effective Radiative Forcing|CO2": "CO2",
    "Effective Radiative Forcing|CH4": "CH4",
    "Effective Radiative Forcing|N2O": "N2O",
    "Effective Radiative Forcing|Stratospheric Water Vapor": "STRAT_H2O",
    "Effective Radiative Forcing|Stratospheric Ozone": "STRAT_O3",
    "Effective Radiative Forcing|Tropospheric Ozone": "TROP_O3",
    "Emissions|CO2": "CO2",
    "Emissions|CH4": "CH4",
    "Emissions|N2O": "N2O",
    # Heat uptake
    "Heat Uptake": "RIB_glob",
    "Heat Content|Ocean": "OHCTOT",
    "Anomalous Radiation": "anomalous_radiation",
    # concentrations
    "Atmospheric Concentrations|CO2": "CO2",
    "Atmospheric Concentrations|CH4": "CH4",
    "Atmospheric Concentrations|N2O": "N2O",
}
forc_sums = ["Aerosols", "Aerosols|Direct Effect"]
fgas_list = [
    "CFC-11",
    "CFC-12",
    "CFC-113",
    "CFC-114",
    "CFC-115",
    "CH3Br",
    "CCl4",
    "CH3CCl3",
    "HCFC-22",
    "HCFC-141b",
    "HCFC-123",
    "HCFC-142b",
    "H-1211",
    "H-1301",
    "H-2402",
    "HFC125",
    "HFC134a",
    "HFC143a",
    "HFC227ea",
    "HFC23",
    "HFC245fa",
    "HFC32",
    "HFC4310mee",
    "C2F6",
    "C6F14",
    "CF4",
    "SF6",
]
ghg_not_fgas = ["CO2", "CH4", "N2O", "TROP_O3", "STRAT_O3", "STRAT_H2O"]
carbon_cycle_outputs = [
    "Biosphere carbon flux",
    "Ocean carbon flux",
    "Airborne fraction CO2",
    "Biosphere carbon pool",
    "Ocean carbon pool",
]


def get_data_from_forc_common(df_temp, variable, v_dict, volc=0, sun=0):
    """
    Get or calculate forcing when dataframe with forcers
    variable and variable dictionary
    If calculating with volcanic and solar forcing
    methods to obtain these need to be supplied
    """
    years = df_temp.Year[:]
    if variable in forc_sums:
        timeseries = np.zeros(len(years))
        for comp, value in v_dict.items():
            if variable in comp and value not in forc_sums:
                timeseries = timeseries + df_temp[value].to_numpy()
    elif variable in ("Fgas", "GHG"):
        timeseries = np.zeros(len(years))
        for comp in fgas_list:
            timeseries = timeseries + df_temp[comp].to_numpy()
        if variable == "GHG":
            for comp in ghg_not_fgas:
                timeseries = timeseries + df_temp[comp].to_numpy()
    elif variable == "Total_forcing+sunvolc":
        timeseries = df_temp["Total_forcing"].to_numpy()
        timeseries = timeseries + volc
        timeseries = timeseries + sun
    else:
        timeseries = df_temp[variable].to_numpy()
    return years, timeseries


def get_data_from_conc(results, variable):
    """
    Get data from concentration files
    """
    df_temp = results["concentrations"]
    years = df_temp.Year[:]
    timeseries = df_temp[variable].to_numpy()  # pylint:disable=unsubscriptable-object
    return years, timeseries


def get_data_from_em(results, variable):
    """
    Get data from emissions files
    """
    df_temp = results["emissions"]
    years = df_temp.Year[:]
    # If concentrations run, CO2 emissions should be taken from carbon cycle back calculation
    if "Emissions" in results["carbon cycle"].keys() and variable == "CO2":
        timeseries = results["carbon cycle"][
            "Emissions"
        ].to_numpy()  # pylint:disable=unsubscriptable-object
    else:
        timeseries = df_temp[
            variable
        ].to_numpy()  # pylint:disable=unsubscriptable-object
    return years, timeseries


def get_data_from_temp_or_rib(results, variable):
    """
    Get data for temperature or rib variables
    """
    return results[variable]


def get_data_from_ohc(results, variable):
    """
    Get data from ocean heat content
    """
    df_temp = results[variable]
    # Units are 10^22J and output should be 10^21J = ZJ
    conv_factor = 10.0
    timeseries = df_temp * conv_factor  # pylint:disable=unsubscriptable-object
    return timeseries


def get_carbon_cycle_outputs(results, variable):
    """
    Get data from carbon cycle output
    """
    if variable.endswith("flux"):
        timeseries = results[variable]
        unit = "Pg C / yr"
    elif variable.endswith("pool"):
        timeseries = np.cumsum(results[variable.replace("pool", "flux")].values)
        unit = "Pg C"
    else:
        timeseries = results[variable]
        unit = "Unitless"
    return timeseries, unit


class CSCMREADER:
    """
    Class to read CICERO-SCM output data
    """

    def __init__(self, nystart, nyend):
        self.variable_dict = openscm_to_cscm_dict
        self.variable_dict["Effective Radiative Forcing|Aerosols|Direct Effect|SOx"] = (
            "SO4_DIR"
        )
        self.temp_list = (
            "dT_glob",
            "dT_glob_air",
            "dT_glob_sea",
            "RIB_glob",
            "anomalous_radiation",
        )
        self.ohc_list = "OHCTOT"
        self.indices = np.arange(nystart, nyend + 1)

    def get_variable_timeseries(self, results, variable, sfilewriter):
        """
        Get variable timeseries
        Connecting up to correct data dictionary to get data
        """
        years, timeseries, unit = (
            pd.Series([], dtype="float64"),
            pd.Series([], dtype="float64"),
            "NoUnit",
        )
        if "Concentration" in variable:
            years, timeseries = get_data_from_conc(
                results, self.variable_dict[variable]
            )
            unit = sfilewriter.concunits[
                sfilewriter.components.index(self.variable_dict[variable])
            ]

        elif "Emissions" in variable:
            years, timeseries = get_data_from_em(results, self.variable_dict[variable])
            unit = sfilewriter.units[
                sfilewriter.components.index(self.variable_dict[variable])
            ]

        elif "Forcing" in variable:
            years, timeseries = self.get_data_from_forc(
                results, self.variable_dict[variable]
            )
            unit = "W/m^2"

        elif variable in carbon_cycle_outputs:
            if "carbon cycle" in results.keys():
                years = self.indices
                timeseries, unit = get_carbon_cycle_outputs(
                    results["carbon cycle"], variable
                )

        elif self.variable_dict[variable] in self.temp_list:
            timeseries = get_data_from_temp_or_rib(
                results, self.variable_dict[variable]
            )
            years = self.indices
            if self.variable_dict[variable] == "RIB_glob":
                unit = "W/m^2"
            else:
                unit = "K"

        elif self.variable_dict[variable] in self.ohc_list:
            timeseries = get_data_from_ohc(results, self.variable_dict[variable])
            years = self.indices
            unit = "ZJ"

        return years, timeseries, unit

    def get_volc_forcing(self, results):
        """
        Return volcanic forcing time series from startyear up to and including endyear
        """
        volc_series = pd.Series(
            (results["Volcanic_forcing_NH"] + results["Volcanic_forcing_SH"]) / 2
        )
        volc_series.index = self.indices  # TODO get correct time rang
        return volc_series

    def get_sun_forcing(self, results):
        """
        Return volcanic forcing time series from startyear up to and including endyear
        """
        sun_series = pd.Series(results["Solar_forcing"])
        sun_series.index = self.indices
        return sun_series

    def get_data_from_forc(self, results, variable):
        """
        Get data from forcing files
        """
        df_temp = results["forcing"]
        if variable == "Total_forcing+sunvolc":
            volc = self.get_volc_forcing(results)
            sun = self.get_sun_forcing(results)
            return get_data_from_forc_common(
                df_temp, variable, self.variable_dict, volc, sun
            )
        years, timeseries = get_data_from_forc_common(
            df_temp, variable, self.variable_dict
        )
        return years, timeseries
