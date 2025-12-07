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
        "efficacy": 1,  # Efficacy of deep ocean heat uptake
        "foan": 0.61,  # Northern hemisphere ocean area fraction
        "foas": 0.81,  # Southern hemisphere ocean area fraction
    }

    output_dict_default = {
        "RIB_glob": "RIB",
        "RIB_N": "RIBN",
        "RIB_S": "RIBS",
        "dT_glob": "dtemp",
        "dT_NH": "dtempnh",
        "dT_SH": "dtempsh",
        "dT_glob_air": "dtemp_air",
        "dT_NH_air": "dtempnh_air",
        "dT_SH_air": "dtempsh_air",
        "dT_glob_sea": "dtemp_sea",
        "dT_NH_sea": "dtempnh_sea",
        "dT_SHsea": "dtempsh_sea",
        "OHC700": "OHC700",
        "OHCTOT": "OHCTOT",
        "dT_fast": "dtemp_fast",
        "dT_slow": "dtemp_slow",
        "OHC_MIXED": "OHC_MIXED",
        "OHC_DEEP": "OHC_DEEP",
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
            - efficacy: Efficacy of deep ocean heat uptake
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
            * self.pamset["efficacy"]
            * (self.temp_fast - self.temp_slow)
        ) / self.pamset["c_fast"]

        # Slow layer temperature change
        dtemp_slow = (
            self.pamset["k"] * (self.temp_fast - self.temp_slow) / self.pamset["c_slow"]
        )

        # Update temperatures
        self.temp_fast += dtemp_fast
        self.temp_slow += dtemp_slow

        # Calculate air temperature response (pure atmospheric response
        # without ocean coupling)
        # Air responds directly to forcing with climate feedback
        forc_nh_air = forc_nh + np.mean(fn_volc)
        forc_sh_air = forc_sh + np.mean(fs_volc)

        # Air temperature = forcing / lambda (simplified energy balance
        # for atmosphere)
        tempn_air = forc_nh_air / self.pamset["lambda"]
        temps_air = forc_sh_air / self.pamset["lambda"]

        # Ocean surface temperature is the fast layer temperature (mixed layer)
        tempn_sea = self.temp_fast
        temps_sea = (
            self.temp_fast
        )  # Two-layer model is globally averaged, so N=S for ocean

        # Combined temperature following upwelling diffusion model pattern:
        # Total = ocean_fraction * ocean_temp + (1 - ocean_fraction) * air_temp
        tempn = (
            self.pamset["foan"] * tempn_sea + (1.0 - self.pamset["foan"]) * tempn_air
        )
        temps = (
            self.pamset["foas"] * temps_sea + (1.0 - self.pamset["foas"]) * temps_air
        )

        # Calculate hemisphere-specific radiative imbalances
        # Following upwelling diffusion model: RIB = forcing - lambda *
        # temperature
        ribn = forc_nh_air - self.pamset["lambda"] * tempn
        ribs = forc_sh_air - self.pamset["lambda"] * temps

        rib_toa = (
            forc
            - self.pamset["lambda"] * self.temp_fast
            - (self.pamset["efficacy"] - 1)
            * self.pamset["k"]
            * (self.temp_fast - self.temp_slow)
        )

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
            "dtemp": (tempn + temps)
            / 2.0,  # Global mean temperature change (weighted combination)
            "dtemp_fast": self.temp_fast,  # Fast layer temperature
            "dtemp_slow": self.temp_slow,  # Slow layer temperature
            "RIB": rib_toa,  # Radiative imbalance at top of atmosphere
            "OHCTOT": ohc_total,  # Total ocean heat content (10^22 J)
            "OHC_MIXED": ohc_mixed,  # Mixed layer ocean heat content (10^22 J)
            "OHC_DEEP": ohc_deep,  # Deep layer ocean heat content (10^22 J)
            # Air/sea temperature components following upwelling diffusion
            # model pattern
            "dtempnh": tempn,  # Northern hemisphere total temperature
            "dtempsh": temps,  # Southern hemisphere total temperature
            "dtemp_air": (tempn_air + temps_air) / 2.0,  # Global air temperature
            "dtempnh_air": tempn_air,  # Northern hemisphere air temperature
            "dtempsh_air": temps_air,  # Southern hemisphere air temperature
            "dtemp_sea": (tempn_sea + temps_sea)
            / 2.0,  # Global sea surface temperature
            "dtempnh_sea": tempn_sea,  # Northern hemisphere sea surface
            # temperature
            "dtempsh_sea": temps_sea,  # Southern hemisphere sea surface
            # temperature
            "RIBN": ribn,  # Northern hemisphere radiative imbalance
            "RIBS": ribs,  # Southern hemisphere radiative imbalance
            "OHC700": ohc700,  # Ocean heat content down to 700m (10^22 J)
        }
