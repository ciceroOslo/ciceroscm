"""
CICEROSCM 
"""
import logging

from .upwelling_diffusion_model import UpwellingDiffusionModel


class CICEROSCM:
    """
    Main ciceroscm class
    """

    def __init__(self):
        """
        Intialise CICEROSCM 
        """
        self.nystart = 1750
        self.nyend = 2100
        self.idtm = 24
        self.scenstart = 1991
        self.emstart = emstart
        self.scenend = self.nyend
        self.emend = self.nyend
        self.lamb = 0.8
        self.qbmb = 0.03
        self.qo3 = 0.34
        self.qdirso2 = -0.4
        self.qindso2 = -0.8
        self.qbc = 0.22
        self.qoc = -0.05

    def _run(self, pamset, cfg):
        """
        Run CICEROSCM
        """
        pass

    def ocean_temperture(TS, TN, dz, sea_fracSH=0.81, sea_fracNH=0.61):
        area_hemi = 2.55e14
        rho = 1030.0
        c = 3997 * 1.0e-22

        oceantemp = (TN * dz * areaNH + TS * dz * areaSH) * rho * c
