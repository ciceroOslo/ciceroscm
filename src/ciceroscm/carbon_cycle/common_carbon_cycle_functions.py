"""
Common carbon cycle functionality
"""

import numpy as np

# Area of the ocean (m^2)
OCEAN_AREA = 3.62e14

# Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
GE_COEFF = 1.0 / (OCEAN_AREA * 9.06)


# Conversion factor ppm/kg --> umol*/m3
PPMKG_TO_UMOL_PER_VOL = 1.722e17

# Conversion factor ppm CO2 -> kg
PPM_CO2_TO_PG_C = 2.123


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
