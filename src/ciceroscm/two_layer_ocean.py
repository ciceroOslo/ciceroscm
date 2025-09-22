"""
2 layer model
"""

import numpy as np


class TwoLayerOceanModel:  # pylint: disable=too-few-public-methods
    """
    Two layer Model with 2 thermal timescales.
    """

    def __init__(self, params = None):
        """
        Initialize with parameters for multiple thermal timescales.

        Parameters
        ----------
        params : dict
            Dictionary containing model parameters:
            - lambda: Climate feedback parameter (W/m^2/K)
            - c_fast: Heat capacity of the fast-response layer (J/m^2/K)
            - c_slow: Heat capacity of the slow-response layer (J/m^2/K)
            - k: Coupling coefficient between fast and slow layers (W/m^2/K)
        """
        if params is None:
            params = {}
        self.lambda_ = params.get("lambda", 3.74 / 3)
        self.c_fast = params.get("mixed", 50) * 1000 * 4181 / (365 * 24 * 3600)
        self.c_slow = (
            params.get("deep", 1200) * 1000 * 4181 / (365 * 24 * 3600)
        )  # Default value: 100.0
        self.k = params.get("k", 0.5)  # eta, Default value: 0.5
        self.efficacy = params.get("efficacy", 1)

        # Initialize temperatures for fast and slow layers
        self.temp_fast = 0.0
        self.temp_slow = 0.0

    def energy_budget(self, forc_nh, forc_sh, fn_volc, fs_volc):
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
            Dictionary containing temperature changes for fast and slow layers.
        """
        forc = (forc_nh + forc_sh) / 2 + np.mean(fn_volc + fs_volc)

        # Fast layer temperature change
        dtemp_fast = (
            forc
            - self.temp_fast * self.lambda_
            - self.k * self.efficacy * (self.temp_fast - self.temp_slow)
        ) / self.c_fast

        # Slow layer temperature change
        dtemp_slow = self.k * (self.temp_fast - self.temp_slow) / self.c_slow

        # Update temperatures
        self.temp_fast += dtemp_fast
        self.temp_slow += dtemp_slow
        rib_toa = (
            forc
            - self.lambda_ * self.temp_fast
            - (self.efficacy - 1) * self.k * (self.temp_fast - self.temp_slow)
        )

        # TODO: Add calculations of Ocean heat content and RIB? Should be knowable/calculable even from this, right?
        return {
            "dtemp": self.temp_fast,  # Global mean temperature change
            "dtemp_fast": self.temp_fast,
            "dtemp_slow": self.temp_slow,
            "dtemp_air": 0.0,
            "dtempnh": 0.0,
            "dtempsh": 0.0,
            "dtempnh_air": 0.0,
            "dtempsh_air": 0.0,
            "dtemp_sea": 0.0,
            "dtempnh_sea": 0.0,
            "dtempsh_sea": 0.0,
            "RIBN": 0.0,
            "RIBS": 0.0,
            "RIB": rib_toa,
            "OHC700": 0.0,
            "OHCTOT": 0.0,
        }
