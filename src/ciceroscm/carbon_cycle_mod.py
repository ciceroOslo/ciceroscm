"""
Stub version of the CarbonCycleModel.
"""
from ._utils import cut_and_check_pamset
import numpy as np
PPM_CO2_TO_PG_C = 2.123

GE_COEFF = 4000

# Conversion factor ppm/kg --> umol*/m3
PPMKG_TO_UMOL_PER_VOL = 1.722e17

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


def take_out_missing(pamset):
    """
    Take out values that are missing from pamset

    Needed to take care of rs_function and rb_function

    Parameters
    ----------
    pamset : dict
        parameter set dictionary which might contain the value "missing"
        for the rs_function and rb_function because of the way we deal
        with expected values in the parameter set. In this case
        we need to take them out of the parameterset so they can be run
        with defaults

    Returns
    -------
    dict
        Updated version of the parameterset with "missing" values deleted
    """
    for key, value in pamset.items():
        if value == "missing":
            del pamset[key]
    return pamset

class CarbonCycleModel:
    """
    Simplified CarbonCycleModel that returns dummy values.
    """

    def __init__(self, pamset):
        """
        Initialise Carbon cycle model

        Parameters
        ----------
            pamset : dict
        """
        pamset = cut_and_check_pamset(
            {
                "idtm": 24,
                "nystart": 1750,
                "nyend": 2100,
                "mixed_carbon": 75.0,
            },
            pamset,
            used={"rs_function": "missing", "rb_function": "missing"},
        )
        self.pamset = take_out_missing(pamset.copy())
        self.pamset["years_tot"] = pamset["nyend"] - pamset["nystart"] + 1
        self.reset_co2_hold()

        
    def reset_co2_hold(self, beta_f=0.287, mixed_carbon=75.0,fnpp_temp_coeff=np.NaN):
        """
        Reset values of CO2_hold for new run

        This method is mainly called to do a new run with the same cscm instance,
        in which case you need to reset hold values, and be able to update
        parameter values for the carbon cycle free parameter beta_f
        """
        self.co2_hold = {
            "yCO2": 278.0,
            "xCO2": 278.0,
            "sCO2": np.zeros(self.pamset["idtm"] * self.pamset["years_tot"]),
            "emCO2_prev": 0.0,
            "dfnpp": np.zeros(self.pamset["idtm"] * self.pamset["years_tot"]),
            "ss1": 0.0,
            "sums": 0.0,
        }
   
    def co2em2conc(self, yr, em_co2_common, dtemp=0):
        """
        Calculate co2 concentrations from emissions, single slab ocean - no land

        Method to calculate co2 concentrations from emissions
        Implementing a rudimentary carbon cycle which loops over
        idtm (usually 24) timesteps a year

        Parameters
        ----------
        yr : int
          Year for which to calculate
        em_co2_common : float
             Sum of CO2 emissions from fossil fuels, land use change and natural emissions
             for the year in question

        Returns
        -------
        float
             CO2 concetrations for year in question
        """
        # TIMESTEP (YR)
        dt = 1.0 / self.pamset["idtm"]

        yr_ix = yr - self.pamset["nystart"]
        # Monthloop:
        for i in range(self.pamset["idtm"]):
            it = yr_ix * self.pamset["idtm"] + i
            sumf = 0.0

  
            # Sum anthropogenic and biospheric and PPMKG_TO_UMOL_PER_VOLert gC/yr --> ppm/yr
            em_co2 = (em_co2_common) / PPM_CO2_TO_PG_C/self.pamset["idtm"]

            self.co2_hold["yco2"] = self.co2_hold["yCO2"]+(self.co2_hold["xCO2"]-self.co2_hold["yCO2"])/GE_COEFF

     
            self.co2_hold["xCO2"] = (
                self.co2_hold["xCO2"] ) + em_co2-(self.co2_hold["xCO2"]-self.co2_hold["yCO2"])/GE_COEFF
            #)
            # print("it: %d, emCO2: %e, sCO2: %e, zCO2: %e, yCO2: %e, xCO2: %e, ss1: %e, ss2: %e, dnfpp:%e"%(it, em_co2, self.co2_hold["sCO2"][it], z_co2, self.co2_hold["yCO2"], self.co2_hold["xCO2"], self.co2_hold["ss1"], ss2, self.co2_hold["dfnpp"][it]))
        return self.co2_hold["xCO2"]

    def get_biosphere_carbon_flux(self, conc_run=False, co2_conc_series=None):
        """
        Return a dummy biosphere carbon flux timeseries.
        """
        return np.zeros(self.pamset["years_tot"])

    def get_ocean_carbon_flux(self, conc_run=False, co2_conc_series=None):
        """
        Return a dummy ocean carbon flux timeseries.
        """
        return np.zeros(self.pamset["years_tot"])

    def back_calculate_emissions(self, co2_conc_series):
        """
        Return a dummy emissions timeseries.
        """
        return np.zeros(len(co2_conc_series))