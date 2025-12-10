"""
2 layer model
"""

import numpy as np

from ..constants import (
    DAY_YEAR,
    OCEAN_AREA,
    SEC_DAY,
    WATER_DENSITY,
    WATER_HEAT_CAPACITY,
)
from .abstract_thermal_model import AbstractThermalModel


class TwoLayerOceanModel(
    AbstractThermalModel
):  # pylint: disable=too-few-public-methods
    """
    Two layer Model with 2 thermal timescales.
    """

    thermal_model_required_pamset = {
        "lambda": 3.74 / 3,  # Climate feedback parameter (W/m^2/K)
        "mixed": 50,  # Ocean mixed layer depth (m)
        "deep": 1200,  # Deep ocean layer depth (m)
        "k": 0.5,  # Coupling coefficient between layers (W/m^2/K)
        "ocean_efficacy": 1,  # Efficacy of deep ocean heat uptake
        "foan": 0.61,  # Northern hemisphere ocean area fraction
        "foas": 0.81,  # Southern hemisphere ocean area fraction
    }

    output_dict_default = {
        "RIB_glob": "RIB",
        "dT_glob": "dtemp",
        "OHC_MIXED": "OHC_MIXED",
        "OHC_DEEP": "OHC_DEEP",
        "OHC700": "OHC700",
        "OHCTOT": "OHCTOT",
        "dT_fast": "dtemp_fast",
        "dT_slow": "dtemp_slow",
    }

    def __init__(self, pamset=None):
        """
        Initialize with parameters for multiple thermal timescales.

        Parameters
        ----------
        pamset : dict
            Dictionary containing model parameters:
            - lambda: Climate feedback parameter (W/m^2/K)
            - mixed: Mixed layer depth (m)
            - deep: Deep ocean layer depth (m)
            - k: Coupling coefficient between layers (W/m^2/K)
            - ocean_efficacy: Efficacy of deep ocean heat uptake
            - foan: Northern hemisphere ocean area fraction
            - foas: Southern hemisphere ocean area fraction
        """
        # Call parent constructor which handles parameter validation and
        # filtering
        super().__init__(pamset)

        # Calculate derived parameters and store in pamset
        # Heat capacity per unit area (J/m^2/K) = depth(m) * density(kg/m^3) *
        # specific_heat(J/kg/K) / time_conversion
        self.pamset["c_fast"] = (
            self.pamset["mixed"]
            * WATER_DENSITY
            * WATER_HEAT_CAPACITY
            / (SEC_DAY * DAY_YEAR)
        )
        self.pamset["c_slow"] = (
            self.pamset["deep"]
            * WATER_DENSITY
            * WATER_HEAT_CAPACITY
            / (SEC_DAY * DAY_YEAR)
        )

        # Initialize temperatures for fast and slow layers
        self.temp_fast = 0.0
        self.temp_slow = 0.0

    def energy_budget(
        self, forc_nh, forc_sh, fn_volc, fs_volc
    ):  # pylint: disable=too-many-locals
        """
        Calculate temperature response with multiple thermal timescales.

        Parameters
        ----------
        forc_nh : float
               Northern hemispheric forcing (W/m^2)
        forc_sh : float
               Southern hemispheric forcing (W/m^2)
        fn_volc : float
               Northern hemispheric volcanic forcing (W/m^2)
        fs_volc : float
               Northern hemispheric volcanic forcing (W/m^2)

        Returns
        -------
        dict
            Dictionary containing temperature changes and diagnostics:
            - dtemp: Global mean temperature (weighted air/sea combination)
            - dtemp_fast: Fast layer (mixed layer) temperature
            - dtemp_slow: Slow layer (deep ocean) temperature
            - dtempnh/dtempsh: Hemispheric total temperatures
            - dtemp_air/dtempnh_air/dtempsh_air: Pure atmospheric temperatures
            - dtemp_sea/dtempnh_sea/dtempsh_sea: Ocean surface temperatures
            - RIB/RIBN/RIBS: Radiative imbalances (global and
              hemisphere-specific)
            - OHCTOT/OHC_MIXED/OHC_DEEP/OHC700: Ocean heat content diagnostics
        """
        forc = (forc_nh + forc_sh) / 2 + np.mean(fn_volc + fs_volc)

        # Fast layer temperature change
        dtemp_fast = (
            forc
            - self.temp_fast * self.pamset["lambda"]
            - self.pamset["k"]
            * self.pamset["ocean_efficacy"]
            * (self.temp_fast - self.temp_slow)
        ) / self.pamset["c_fast"]

        # Slow layer temperature change
        dtemp_slow = (
            self.pamset["k"] * (self.temp_fast - self.temp_slow) / self.pamset["c_slow"]
        )

        # Update temperatures
        self.temp_fast += dtemp_fast
        self.temp_slow += dtemp_slow

        # TODO: Clean up and just calculate dtemp

        rib_toa = (
            forc
            - self.pamset["lambda"] * self.temp_fast
            - (self.pamset["ocean_efficacy"] - 1)
            * self.pamset["k"]
            * (self.temp_fast - self.temp_slow)
        )
        dtemp = (
            forc
            / self.pamset["lambda"]
            * (1 - (self.pamset["foan"] + self.pamset["foas"]) / 2.0)
            + self.temp_fast * (self.pamset["foan"] + self.pamset["foas"]) / 2
        )  # Approximate global mean temperature change

        # Calculate ocean heat content (10^22 J)
        # OHC = temperature_change * heat_capacity_per_unit_area *
        # ocean_area / 1e22
        # Following UDM convention: convert from J/m² to 10^22 J units

        ohc_mixed = (
            self.temp_fast
            * self.pamset["c_fast"]
            * (SEC_DAY * DAY_YEAR)
            * OCEAN_AREA
            / 1e22
        )  # Convert to 10^22 J
        ohc_deep = (
            self.temp_slow
            * self.pamset["c_slow"]
            * (SEC_DAY * DAY_YEAR)
            * OCEAN_AREA
            / 1e22
        )  # Convert to 10^22 J
        ohc_total = ohc_mixed + ohc_deep

        # Calculate OHC700 (ocean heat content down to 700m)
        # For 2-layer model, we need to estimate how much of each layer
        # contributes to the top 700m
        mixed_depth = self.pamset["mixed"]  # Mixed layer depth (m)
        deep_depth = self.pamset["deep"]  # Deep layer depth (m)

        if mixed_depth >= 700:
            # If mixed layer is deeper than 700m, OHC700 is just a
            # fraction of mixed layer
            ohc700 = ohc_mixed * (700 / mixed_depth)
        else:
            # Mixed layer + part of deep layer contributes to top 700m
            remaining_depth = 700 - mixed_depth
            deep_fraction = min(remaining_depth / deep_depth, 1.0)
            ohc700 = ohc_mixed + ohc_deep * deep_fraction

        # Return outputs: meaningful values for this model structure plus
        # compatibility placeholders for interface consistency
        return {
            "dtemp": dtemp,  # Global mean temperature change
            "dtemp_fast": self.temp_fast,  # Fast layer temperature
            "dtemp_slow": self.temp_slow,  # Slow layer temperature
            "RIB": rib_toa,  # Radiative imbalance at top of atmosphere
            "OHCTOT": ohc_total,  # Total ocean heat content (10^22 J)
            "OHC_MIXED": ohc_mixed,  # Mixed layer ocean heat content (10^22 J)
            "OHC_DEEP": ohc_deep,  # Deep layer ocean heat content (10^22 J)
            "dtemp_sea": self.temp_fast,  # Global sea surface temperature
            "OHC700": ohc700,  # Ocean heat content down to 700m (10^22 J)
        }
