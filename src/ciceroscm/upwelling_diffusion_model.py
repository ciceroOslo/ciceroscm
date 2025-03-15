"""
Stub version of the UpwellingDiffusionModel.
"""

import numpy as np

class UpwellingDiffusionModel:
    """
    Simplified UpwellingDiffusionModel that returns dummy values.
    """


    def __init__(self, params):
        """
        Intialise

        Setting up heights first, then starting up empty dict and
        arrays , and some starting values

        Parameters
        ----------
        params : dict
              Physical parameters to define the instance
        """

        # Dict of values to keep from one year to the next
        self.prev_values = {
            "dtemp": 0.0,
        }

        self.dtempprev = 0.0

    def energy_budget(self, forc_nh, forc_sh, fn_volc, fs_volc):
        """
        Return dummy energy budget results.
        """
        forc=forc_nh+forc_sh
        dtemp = self.prev_values['dtemp']#+(forc-self.prev_values['dtemp']/1.1)/1e20
        self.prev_values['dtemp'] = dtemp
        return {
            "dtemp": dtemp,
            "dtempnh": 0.0,
            "dtempsh": 0.0,
            "dtemp_air": 0.0,
            "dtempnh_air": 0.0,
            "dtempsh_air": 0.0,
            "dtemp_sea": 0.0,
            "dtempnh_sea": 0.0,
            "dtempsh_sea": 0.0,
            "RIBN": 0.0,
            "RIBS": 0.0,
            "RIB": 0.0,
            "OHC700": 0.0,
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