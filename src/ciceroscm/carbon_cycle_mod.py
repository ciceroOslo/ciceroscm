"""
Simple CarbonCycleModel with temperature feedback and no pre-computed pulse response functions.
"""

import numpy as np

PPM_CO2_TO_PG_C = 2.123  # Conversion factor ppm CO2 -> PgC
OCEAN_AREA = 3.62e14  # Area of the ocean in m^2
GE_COEFF = 1.0 / (OCEAN_AREA * 9.06)  # Gas exchange coefficient (yr^-1*m^-2)


def calculate_airborne_fraction(em_timeseries, conc_timeseries):
    """
    Calculate Airborne Fraction of CO2 from emissions timeseries

    Parameters
    ----------
    em_timeseries: np.ndarray
        Emissions timeseries, either inputs or backcalculated for concentration
        run
    conc_timeseries : np.ndarray
        Concentrations timeseries. Should be the same length as the emissions
        timeseries

    Returns
    -------
    np.ndarray
        Airborne fraction calculated from the em_timeseries and conc_timeseries
    """
    airborne_fraction = (
        (conc_timeseries - 278.0) / np.cumsum(em_timeseries) * PPM_CO2_TO_PG_C
    )
    return airborne_fraction


class CarbonCycleModel:
    """
    CarbonCycleModel with explicit carbon pools and temperature-dependent dynamics.
    """

    def __init__(self, pamset):
        """
        Initialize the Carbon Cycle Model.

        Parameters
        ----------
        pamset : dict
            Parameter set for the model, containing:
            - idtm: Number of subyearly timesteps (e.g., 24 for monthly steps).
            - nystart: Start year of the simulation.
            - nyend: End year of the simulation.
            - beta_f: CO2 fertilization factor (affects land carbon uptake).
            - land_temp_sensitivity: Sensitivity of NPP to temperature (PgC/K).
            - soil_respiration_rate: Rate of soil carbon decay (yr^-1).
            - ocean_mixed_layer_depth: Depth of the ocean mixed layer (m).
            - ocean_exchange_rate: Rate of carbon exchange between mixed layer and deep ocean (yr^-1).
            - vegetation_to_soil_fraction: Fraction of vegetation carbon transferred to soil per year.
            - ocean_solubility_base: Base solubility of CO2 in the ocean (PgC/ppm).
            - ocean_solubility_temp_coeff: Temperature sensitivity of ocean CO2 solubility.
        """
        self.pamset = {
            "idtm": pamset.get("idtm", 24),
            "nystart": pamset.get("nystart", 1750),
            "nyend": pamset.get("nyend", 2100),
            "beta_f": pamset.get("beta_f", 0.287),
            "land_temp_sensitivity": pamset.get("land_temp_sensitivity", .1),
            "soil_respiration_rate": pamset.get("soil_respiration_rate", 0.02),
            "ocean_mixed_layer_depth": pamset.get("ocean_mixed_layer_depth", 25.0),
            "ocean_exchange_rate": pamset.get("ocean_exchange_rate", 0.01),
            "vegetation_to_soil_fraction": pamset.get("vegetation_to_soil_fraction", 0.1),
            "ocean_solubility_base": pamset.get("ocean_solubility_base", 0.02),
            "ocean_solubility_temp_coeff": pamset.get("ocean_solubility_temp_coeff", -0.01),
        }
        self.pamset["years_tot"] = self.pamset["nyend"] - self.pamset["nystart"] + 1

        # Initialize carbon pools
        self.atmospheric_co2 = 278.0  # Pre-industrial CO2 concentration (ppm)
        self.vegetation_carbon = 600.0  # Vegetation carbon pool (PgC)
        self.soil_carbon = 3000.0  # Soil carbon pool (PgC)
        self.ocean_mixed_layer_carbon = 0.0  # Ocean mixed layer carbon anomaly (PgC)
        self.ocean_deep_carbon = 0.0  # Deep ocean carbon anomaly (PgC)

    def co2em2conc(self, yr, em_co2_common, dtemp=0):
        """
        Calculate CO2 concentrations from emissions.

        Parameters
        ----------
        yr : int
            Year for which to calculate.
        em_co2_common : float
            Total CO2 emissions (GtC/yr).
        dtemp : float
            Global mean temperature difference from pre-industrial (K).

        Returns
        -------
        float
            Updated atmospheric CO2 concentration (ppm).
        """
        
        dt = 1.0# / self.pamset["idtm"]  # Timestep length (years)

        # Convert emissions from GtC/yr to PgC/yr
        emissions = em_co2_common

        # Land carbon fluxes
        npp = self.calculate_npp(dtemp)  # Net Primary Production (PgC/yr)
        vegetation_to_soil = self.vegetation_carbon * self.pamset["vegetation_to_soil_fraction"]  # Fraction of vegetation carbon transferred to soil
        soil_respiration = self.soil_carbon * self.pamset["soil_respiration_rate"]  # Soil carbon decay (PgC/yr)

        # Update land carbon pools
        self.vegetation_carbon += (npp - vegetation_to_soil) * dt
        self.soil_carbon += (vegetation_to_soil - soil_respiration) * dt

        # Ocean carbon fluxes
        ocean_uptake = self.calculate_ocean_uptake(dtemp)  # Ocean carbon uptake (PgC/yr)
        mixed_to_deep = self.ocean_mixed_layer_carbon * self.pamset["ocean_exchange_rate"]  # Mixed to deep ocean flux
        deep_to_mixed = self.ocean_deep_carbon * self.pamset["ocean_exchange_rate"]  # Deep to mixed ocean flux

        # Update ocean carbon pools
        self.ocean_mixed_layer_carbon += (ocean_uptake - mixed_to_deep + deep_to_mixed) * dt
        self.ocean_deep_carbon += (mixed_to_deep - deep_to_mixed) * dt

        # Net flux to atmosphere
        net_flux = emissions - ocean_uptake + soil_respiration - npp

        # Update atmospheric CO2 concentration (ppm)
        self.atmospheric_co2 += net_flux / PPM_CO2_TO_PG_C * dt
        return self.atmospheric_co2

    def calculate_npp(self, dtemp):
        """
        Calculate Net Primary Production (NPP) as a function of temperature and CO2.

        Parameters
        ----------
        dtemp : float
            Global mean temperature difference from pre-industrial (K).

        Returns
        -------
        float
            NPP (PgC/yr).
        """
        co2_fertilization = self.pamset["beta_f"] * np.log(self.atmospheric_co2 / 278.0)
        temp_effect = self.pamset["land_temp_sensitivity"] * dtemp
        return max(0.0, 60.0 + co2_fertilization - temp_effect)  # Ensure NPP is non-negative

    def calculate_ocean_uptake(self, dtemp):
        """
        Calculate ocean carbon uptake as a function of temperature and CO2.

        Parameters
        ----------
        dtemp : float
            Global mean temperature difference from pre-industrial (K).

        Returns
        -------
        float
            Ocean carbon uptake (PgC/yr).
        """
        solubility = self.pamset["ocean_solubility_base"] * (
            1 + self.pamset["ocean_solubility_temp_coeff"] * dtemp
        )  # Solubility decreases with warming
        return solubility * (self.atmospheric_co2 - 278.0)  # Uptake proportional to CO2 difference
    
    def reset_co2_hold(self, beta_f=0.287, mixed_carbon=75.0, fnpp_temp_coeff=0):
        """
        stub
        """
