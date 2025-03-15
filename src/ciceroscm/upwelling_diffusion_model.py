"""
Simplest possible thermal model
"""

import numpy as np

class UpwellingDiffusionModel:
    """
    Highly Simplified UpwellingDiffusionModel 
    """

    def __init__(self, params):
        """
        Initialize with minimal parameters.
        """
        self.pamset = {
            "lm": 40,
            "ldtime": 12,
        }
        self.dz = np.ones(self.pamset["lm"]) * 100.0
        self.tn = np.zeros(self.pamset["lm"])
        self.ts = np.zeros(self.pamset["lm"])
        self.fdb = 1.0 / params["lambda"]
        self.c1 = 1/500
        self.prev_values = {
            "fn": 0.0,
            "fs": 0.0,
            "dtemp": 0.0,
        }

        self.dtempprev = 0.0

    def energy_budget(self, forc_nh, forc_sh, fn_volc, fs_volc):
        """
        Return energy budget results.
        """
        forc=forc_nh+forc_sh
        dtemp=self.dtempprev+(forc-self.dtempprev/self.fdb)*self.c1
        self.dtempprev=dtemp
        return {
            "dtemp": dtemp,
            "dtempnh": dtemp,
            "dtempsh": dtemp,
            "dtemp_air": 0.0,
            "dtempnh_air": 0.0,
            "dtempsh_air": 0.0,
            "dtemp_sea": 0.0,
            "dtempnh_sea": 0.0,
            "dtempsh_sea": 0.0,
            "RIBN": 0.0,
            "RIBS": 0.0,
            "RIB": 0.0,
            "OHC700": dtemp*self.c1,
            "OHCTOT": 0.0,
        }

    def ocean_temperature(self):
        """
        Return dummy ocean temperature results.
        """
        return {
            "OHC700": 0.0,
            "OHCTOT": 0.0,
        }