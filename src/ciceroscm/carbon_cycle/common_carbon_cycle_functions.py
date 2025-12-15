"""
Common carbon cycle functionality
"""

import numpy as np

from .._utils import cut_and_check_pamset

# Area of the ocean (m^2)
OCEAN_AREA = 3.62e14

# Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
GE_COEFF = 1.0 / (OCEAN_AREA * 9.06)


# Conversion factor ppm/kg --> umol*/m3
PPMKG_TO_UMOL_PER_VOL = 1.722e17

# Conversion factor ppm CO2 -> kg
PPM_CO2_TO_PG_C = 2.123

PREINDUSTRIAL_CO2_CONC = 278.0


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
    if any(np.cumsum(em_timeseries) == 0):
        return np.ones_like(em_timeseries) * np.nan
    airborne_fraction = (
        (conc_timeseries - PREINDUSTRIAL_CO2_CONC)
        / np.cumsum(em_timeseries)
        * PPM_CO2_TO_PG_C
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


def carbon_cycle_init_pamsets(
    pamset_emiconc, pamset_carbon, carbon_cycle_required_pamset, used=None
):
    """
    Check and combine parameter sets at initilisation of carbon cycle modules

    Parameters
    ----------
    pamset_emiconc : dict
        Parameter set from the concentrations emission handler, containing:
        - idtm: Number of subyearly timesteps (e.g., 24 for monthly steps).
        - nystart: Start year of the simulation.
        - nyend: End year of the simulation.
    pamset_carbon : dict
        Optional carbon specific parameter set allowed options
    carbon_cycle_required_pamset : dict
        Required carbon cycle dictionary content, the values of which
        the optional parameters from pamset_carbon can overwrite
    used : dict
        Optional parameters in pamset_carbon which will be cut if not
        present in pamset_carbon
    """
    pamset = cut_and_check_pamset(
        {
            "idtm": 24,
            "nystart": 1750,
            "nyend": 2100,
        },
        pamset_emiconc,
    )
    if pamset_carbon is None:
        pamset_carbon = carbon_cycle_required_pamset
    else:
        pamset_carbon = cut_and_check_pamset(
            carbon_cycle_required_pamset,
            pamset_carbon,
            used=used,
        )
    if used is not None:
        pamset_carbon = take_out_missing(pamset_carbon.copy())

    return pamset, pamset_carbon
