"""
Physical and mathematical constants used throughout CICERO-SCM

This module contains fundamental physical constants, conversion factors,
and commonly used values across the climate model components.

Note: This is a new centralized constants module. Some constants may still
exist in other modules during the transition period.
"""

# Time conversion constants
SEC_DAY = 86400  # Seconds per day
DAY_YEAR = 365.0  # Days per year (simplified, no leap years)

# Water properties (seawater)
WATER_DENSITY = 1000  # kg/m^3 - Density of seawater
WATER_HEAT_CAPACITY = 4184  # J/kg/K - Specific heat capacity of seawater

# Earth surface areas
OCEAN_AREA = 3.61e14  # m^2 - Global ocean surface area

# Derived constants
YEAR_IN_SECONDS = DAY_YEAR * SEC_DAY  # Seconds per year

# Pattern-mediated feedback (Tier 3): tracers grouped as "aerosol" for the
# magnitude-weighted forcing fraction w_aero(t) that modulates the climate
# feedback parameter. See notebooks/variablelambda/ for derivation.
AEROSOL_TRACERS = (
    "SO2",  # stored as SO2 internally; renamed to SO4_DIR in output only
    "SO4_IND",
    "BC",
    "OC",
    "NOx",
    "NMVOC",
    "NH3",
    "BMB_AEROS",
)
