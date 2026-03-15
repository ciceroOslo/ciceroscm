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
    "Surface Air Ocean Temperature Change": "dT_glob_sea",
    # ERFs
    "Effective Radiative Forcing": "Total_forcing+sunvolc",
    "Effective Radiative Forcing|Anthropogenic": "Total_forcing",
    "Effective Radiative Forcing|Anthropogenic|Aerosol": "Aerosols",
    "Effective Radiative Forcing|Anthropogenic|Aerosol|Aerosol-radiation Interactions": "Aerosols|Direct Effect",
    "Effective Radiative Forcing|Anthropogenic|Aerosol|Aerosol-radiation Interactions|BC": "BC",
    "Effective Radiative Forcing|Anthropogenic|Aerosol|Aerosol-radiation Interactions|OC": "OC",
    "Effective Radiative Forcing|Anthropogenic|Aerosol|Aerosol-radiation Interactions|Sulfate": "SO4_DIR",
    "Effective Radiative Forcing|Anthropogenic|Aerosol|Aerosol-cloud Interactions": "SO4_IND",
    "Effective Radiative Forcing|Greenhouse Gases": "GHG",
    "Effective Radiative Forcing|Anthropogenic|F-Gases": "Fgas",
    "Effective Radiative Forcing|Anthropogenic|CO2": "CO2",
    "Effective Radiative Forcing|Anthropogenic|CH4": "CH4",
    "Effective Radiative Forcing|Anthropogenic|N2O": "N2O",
    "Effective Radiative Forcing|Anthropogenic|Other|Stratospheric H2O": "STRAT_H2O",
    "Effective Radiative Forcing|Anthropogenic|Ozone|Stratospheric Contribution": "STRAT_O3",
    "Effective Radiative Forcing|Anthropogenic|Ozone|Tropospheric Contribution": "TROP_O3",
    "Emissions|CO2": "CO2",
    "Emissions|CH4": "CH4",
    "Emissions|N2O": "N2O",
    # Heat uptake
    "Heat Uptake": "RIB_glob",
    "Heat Content|Ocean": "OHCTOT",
    "Heat Content|Ocean|0-700m": "OHC700",
    "Anomalous Radiation": "anomalous_radiation",
    # concentrations
    "Atmospheric Concentrations|CO2": "CO2",
    "Atmospheric Concentrations|CH4": "CH4",
    "Atmospheric Concentrations|N2O": "N2O",
}
forc_sums = ["Aerosols", "Aerosols|Direct Effect"]
fgas_list = [
    "HFC125",
    "HFC134a",
    "HFC143a",
    "HFC152a",
    "HFC227ea",
    "HFC236fa",
    "HFC23",
    "HFC245fa",
    "HFC32",
    "HFC365mfc",
    "HFC4310mee",
    "NF3",
    "C2F6",
    "C3F8",
    "C4F10",
    "C5F12",
    "C6F14",
    "C7F16",
    "C8F18",
    "cC4F8",
    "CF4",
    "SF6",
    "SO2F2",
]
montreal_list = [
    "CFC-11",
    "CFC-12",
    "CFC-113",
    "CFC-114",
    "CFC-115",
    "CH3Br",
    "CCl4",
    "CH3CCl3",
    "H-1202",
    "H-1211",
    "H-1301",
    "H-2402",
    "HCFC-123",
    "HCFC-142b",
    "HCFC-141b",
    "HCFC-22",
    "CHCl3",
    "CH3Cl",
    "CH2Cl2",
]


def fgas_transform(fgas_name):
    """
    Reformat fgas names to match RCMIP output

    Parameters
    ----------
    fgas_name : str
        Name of fgas to reformat

    Returns
    -------
    str
        Reformatted fgas name
    """
    if fgas_name.startswith("HFC"):
        return f"HFC|{fgas_name}"
    if fgas_name.lower().startswith("c"):
        return f"PFC|{fgas_name}"
    return fgas_name


def montreal_transform(montgas):
    """
    Reformat montreal gas names to match RCMIP output

    Parameters
    ----------
    montgas : str
        Name of montreal gas to reformat

    Returns
    -------
    str
        Reformatted fgas name
    """
    if montgas.startswith("CFC"):
        return f"CFC|{montgas.replace('-', '')}"
    if montgas.startswith("H-"):
        return montgas.replace("H-", "Halon")
    if montgas.startswith("HCFC-"):
        return montgas.replace("HCFC-", "HCFC")
    return montgas


for fgas in fgas_list:
    openscm_to_cscm_dict[
        "Effective Radiative Forcing|Anthropogenic|F-Gases|" + fgas_transform(fgas)
    ] = fgas
    openscm_to_cscm_dict[
        "Atmospheric Concentrations|F-Gases|" + fgas_transform(fgas)
    ] = fgas

for mont in montreal_list:
    openscm_to_cscm_dict[
        "Effective Radiative Forcing|Anthropogenic|Montreal Gases|"
        + montreal_transform(mont)
    ] = mont
    openscm_to_cscm_dict[
        "Atmospheric Concentrations|Montreal Gases|" + montreal_transform(mont)
    ] = mont

ghg_not_fgas = ["CO2", "CH4", "N2O", "TROP_O3", "STRAT_O3", "STRAT_H2O"]
carbon_cycle_outputs = {
    "Carbon Flux|Land": "Biosphere carbon flux",
    "Carbon Flux|Ocean": "Ocean carbon flux",
    "Airborne fraction CO2": "Airborne fraction CO2",
    "Carbon Pool|Land": "Biosphere carbon pool",
    "Carbon Pool|Ocean": "Ocean carbon pool",
}


def get_data_from_forc_common(
    df_temp, variable, v_dict, volc=0, sun=0
):  # pylint:disable=too-many-branches
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
            if comp in df_temp.columns:
                timeseries = timeseries + df_temp[comp].to_numpy()
        if variable == "GHG":
            for comp in montreal_list:
                if comp in df_temp.columns:
                    timeseries = timeseries + df_temp[comp].to_numpy()
            for comp in ghg_not_fgas:
                if comp in df_temp.columns:
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
        self.ohc_list = ["OHCTOT", "OHC700"]
        self.indices = np.arange(nystart, nyend + 1)

    def get_variable_timeseries(
        self, results, variable, sfilewriter
    ):  # pylint:disable=too-many-branches
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
            if variable in self.variable_dict:
                years, timeseries = get_data_from_conc(
                    results, self.variable_dict[variable]
                )
                unit = sfilewriter.concunits[
                    sfilewriter.components.index(self.variable_dict[variable])
                ]

        elif "Emissions" in variable:
            if variable in self.variable_dict:
                years, timeseries = get_data_from_em(
                    results, self.variable_dict[variable]
                )
                unit = sfilewriter.units[
                    sfilewriter.components.index(self.variable_dict[variable])
                ]

        elif "Forcing" in variable:
            if variable in self.variable_dict:
                years, timeseries = self.get_data_from_forc(
                    results, self.variable_dict[variable]
                )
                unit = "W/m^2"

        elif variable in carbon_cycle_outputs:
            if "carbon cycle" in results.keys():
                years = self.indices
                timeseries, unit = get_carbon_cycle_outputs(
                    results["carbon cycle"], carbon_cycle_outputs[variable]
                )

        elif (
            variable in self.variable_dict
            and self.variable_dict[variable] in self.temp_list
        ):
            timeseries = get_data_from_temp_or_rib(
                results, self.variable_dict[variable]
            )
            years = self.indices
            if self.variable_dict[variable] == "RIB_glob":
                unit = "W/m^2"
            else:
                unit = "K"

        elif (
            variable in self.variable_dict
            and self.variable_dict[variable] in self.ohc_list
        ):
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
