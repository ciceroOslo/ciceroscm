import numpy as np

class CarbonCycleModel:
    """
    Simplified CarbonCycleModel that returns dummy values.
    """

    def __init__(self, pamset):
        """
        Initialize with a minimal parameter set.
        """
        self.pamset = {
            "idtm": 24,
            "nystart": 1750,
            "nyend": 2100,
            "years_tot": 2100 - 1750 + 1,
        }

    def co2em2conc(self, yr, em_co2_common):
        """
        Return a dummy CO2 concentration value.
        """
        return 285.0

    def back_calculate_emissions(self, co2_conc_series):
        """
        Return a dummy emissions timeseries.
        """
        return np.zeros(len(co2_conc_series))
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

