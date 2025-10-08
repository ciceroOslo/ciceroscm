"""
Simple CarbonCycleModel with temperature feedback and no pre-computed pulse response functions.
"""

import numpy as np
import pandas as pd

from .carbon_cycle_abstract import AbstractCarbonCycleModel
from .common_carbon_cycle_functions import (  # OCEAN_AREA, GE_COEFF
    PPM_CO2_TO_PG_C,
)


class CarbonCycleModel(AbstractCarbonCycleModel):
    """
    CarbonCycleModel with explicit carbon pools and temperature-dependent dynamics.
    """

    carbon_cycle_model_required_pamset = {
        "beta_f": 0.287,
        "land_temp_sensitivity": 0.1,
        "soil_respiration_rate": 0.02,
        "ocean_mixed_layer_depth": 25.0,
        "ocean_exchange_rate": 0.01,
        "vegetation_to_soil_fraction": 0.1,
        "ocean_solubility_base": 0.02,
        "ocean_solubility_temp_coeff": -0.01,
    }

    def __init__(self, pamset_emiconc, pamset_carbon=None):
        """
        Initialize the Carbon Cycle Model.

        Parameters
        ----------
        pamset_emiconc : dict
            Parameter set from the concentrations emission handler, containing:
            - idtm: Number of subyearly timesteps (e.g., 24 for monthly steps).
            - nystart: Start year of the simulation.
            - nyend: End year of the simulation.
        pamset_carbon : dict
            Optional carbon specific parameter set
            - beta_f: CO2 fertilization factor (affects land carbon uptake).
            - land_temp_sensitivity: Sensitivity of NPP to temperature (PgC/K).
            - soil_respiration_rate: Rate of soil carbon decay (yr^-1).
            - ocean_mixed_layer_depth: Depth of the ocean mixed layer (m).
            - ocean_exchange_rate: Rate of carbon exchange between mixed layer and deep ocean (yr^-1).
            - vegetation_to_soil_fraction: Fraction of vegetation carbon transferred to soil per year.
            - ocean_solubility_base: Base solubility of CO2 in the ocean (PgC/ppm).
            - ocean_solubility_temp_coeff: Temperature sensitivity of ocean CO2 solubility.
        """
        super().__init__(pamset_emiconc, pamset_carbon=pamset_carbon)

    def co2em2conc(
        self, yr, em_co2_common, feedback_dict=None
    ):  # pylint: disable=unused-argument
        """
        Calculate CO2 concentrations from emissions.

        Parameters
        ----------
        yr : int
            Year for which to calculate.
        em_co2_common : float
            Total CO2 emissions (GtC/yr).
        feedback_dict : dict
            Dictionary containing feedback variables (e.g., {"dtemp": float}).

        Returns
        -------
        float
            Updated atmospheric CO2 concentration (ppm).
        """
        dt = 1.0  # / self.pamset["idtm"]  # Timestep length (years)
        if feedback_dict is None:
            dtemp = 0.0
        else:
            dtemp = feedback_dict.get("dtemp", 0.0)

        # Convert emissions from GtC/yr to PgC/yr
        emissions = em_co2_common

        # Land carbon fluxes
        npp = self.calculate_npp(dtemp)  # Net Primary Production (PgC/yr)
        vegetation_to_soil = (
            self.vegetation_carbon * self.pamset["vegetation_to_soil_fraction"]
        )  # Fraction of vegetation carbon transferred to soil
        soil_respiration = (
            self.soil_carbon * self.pamset["soil_respiration_rate"]
        )  # Soil carbon decay (PgC/yr)

        # Update land carbon pools
        self.vegetation_carbon += (npp - vegetation_to_soil) * dt
        self.soil_carbon += (vegetation_to_soil - soil_respiration) * dt

        # Ocean carbon fluxes
        ocean_uptake = self.calculate_ocean_uptake(
            dtemp
        )  # Ocean carbon uptake (PgC/yr)
        mixed_to_deep = (
            self.ocean_mixed_layer_carbon * self.pamset["ocean_exchange_rate"]
        )  # Mixed to deep ocean flux
        deep_to_mixed = (
            self.ocean_deep_carbon * self.pamset["ocean_exchange_rate"]
        )  # Deep to mixed ocean flux

        # Update ocean carbon pools
        self.ocean_mixed_layer_carbon += (
            ocean_uptake - mixed_to_deep + deep_to_mixed
        ) * dt
        self.ocean_deep_carbon += (mixed_to_deep - deep_to_mixed) * dt

        # Net flux to atmosphere
        net_flux = emissions - ocean_uptake + soil_respiration - npp

        # Update atmospheric CO2 concentration (ppm)
        self.atmospheric_co2 += net_flux / PPM_CO2_TO_PG_C * dt
        return self.atmospheric_co2

    def calculate_npp(self, dtemp, co2_conc_series=None):
        """
        Calculate Net Primary Production (NPP) as a function of temperature and CO2.

        Parameters
        ----------
        dtemp : float or np.ndarray
            Global mean temperature difference from pre-industrial (K).
        co2_conc_series : np.ndarray
            Timeseries of atmospheric CO2 concentrations over which to estimate
            npp values

        Returns
        -------
        float
            NPP (PgC/yr).
        """
        temp_effect = self.pamset["land_temp_sensitivity"] * dtemp

        if co2_conc_series is None:
            co2_fertilization = self.pamset["beta_f"] * np.log(
                self.atmospheric_co2 / 278.0
            )
            # TODO: Do we need some error handling here?
            return max(
                0.0, 60.0 + co2_fertilization - temp_effect
            )  # Ensure NPP is non-negative

        co2_fertilization = self.pamset["beta_f"] * np.log(co2_conc_series / 278.0)
        return np.clip(60.0 + co2_fertilization - temp_effect, a_min=0, a_max=None)

    def calculate_ocean_uptake(self, dtemp, co2_conc_series=None):
        """
        Calculate ocean carbon uptake as a function of temperature and CO2.

        Parameters
        ----------
        dtemp : float or np.ndarray
            Global mean temperature difference from pre-industrial (K).
        co2_conc_series : np.ndarray
            Time series of concentrations to calculate the ocean uptake
            over full concnentration time series. When this is sent,
            dtemp is assumed to be a corresponding temperature time series

        Returns
        -------
        float
            Ocean carbon uptake (PgC/yr).
        """
        solubility = self.pamset["ocean_solubility_base"] * (
            1 + self.pamset["ocean_solubility_temp_coeff"] * dtemp
        )  # Solubility decreases with warming
        if co2_conc_series is None:
            return solubility * (
                self.atmospheric_co2 - 278.0
            )  # Uptake proportional to CO2 difference
        # TODO: Do we need some error handling here?
        return solubility * (co2_conc_series - 278.0)

    def reset_co2_hold(self, pamset_carbon=None):
        """
        Reset values of co2 pools

        This method is mainly called to do a new run with the same cscm instance,
        in which case you need to reset pool values, and be able to update
        parameter values for the carbon cycle free parameters

        Parameters
        ----------
        pamset_carbon : dict
            Optional dictionary of new values for a subset of the free
            carbon cycle parameters
        """
        self.atmospheric_co2 = 278.0  # Pre-industrial CO2 concentration (ppm)
        self.vegetation_carbon = 600.0  # Vegetation carbon pool (PgC)
        self.soil_carbon = 3000.0  # Soil carbon pool (PgC)
        self.ocean_mixed_layer_carbon = 0.0  # Ocean mixed layer carbon anomaly (PgC)
        self.ocean_deep_carbon = 0.0  # Deep ocean carbon anomaly (PgC)
        super().reset_co2_hold(pamset_carbon)

    # TODO: Improve
    def get_carbon_cycle_output(
        self, years, conc_run=False, conc_series=None, feedback_dict_series=None
    ):  # pylint: disable=unused-argument
        """
        Get carbon cycle output to print out
        """
        # print(conc_series)
        # print(years)
        if conc_series is None:
            return None
        if feedback_dict_series is None:
            dtemp_series = np.zeros_like(conc_series)
        else:
            dtemp_series = feedback_dict_series.get("dtemp", np.zeros_like(conc_series))
        df_carbon = pd.DataFrame(
            data={
                "Net primary production": self.calculate_npp(dtemp_series, conc_series),
                "Ocean carbon flux": self.calculate_ocean_uptake(
                    dtemp_series, conc_series
                ),
            },
            index=years,
        )
        # print(df_carbon)
        return df_carbon

    def back_calculate_emissions(
        self, co2_conc_series, feedback_dict_series=None
    ):  # pylint: disable=unused-argument
        """
        Do backwards calculation to get emissions time series given
        concentrations and temperature series
        """
        return None
