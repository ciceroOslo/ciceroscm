"""
Physical and mathematical constants used throughout CICERO-SCM

This module contains fundamental physical constants, conversion factors,
and commonly used values across the climate model components.

Constants are organized by category and include documentation
about their usage and derivation where applicable.
"""

# Time conversion constants
SECONDS_PER_DAY = 86400  # Seconds per day
DAYS_PER_YEAR = 365.0  # Days per year (simplified, no leap years)

# Water properties (seawater)
WATER_HEAT_CAPACITY = 4181  # J/kg/K - Specific heat capacity of seawater

# Ocean model specific constants (UDM)
RHO_SEAWATER = 1030.0  # kg/m^3 - Seawater density used in both UDM and Two-Layer models
RHO_OHC = 1.03  # g/cm^3 - Density for UDM heat capacity calculations (legacy value)
CP_UDM = 0.955  # Dimensionless heat capacity factor for UDM
CONV_FAC_UDM = 0.485  # Dimensionless conversion factor for UDM heat capacity
C_OHC = 3.997e-19  # Conversion constant for ocean heat content calculations (units: see UDM docs)
OCEAN_AREA_HEMISPHERE = 2.55e14  # m^2 - Ocean area per hemisphere (used in UDM)
DZ = 100.0  # m - Standard ocean layer thickness in UDM
INDEX_700M = 7  # Layer index constant for 700m depth calculations in UDM

# Earth surface areas
OCEAN_AREA = 3.61e14  # m^2 - Global ocean surface area (used in Two-Layer model)

# Derived constants
YEAR_IN_SECONDS = DAYS_PER_YEAR * SECONDS_PER_DAY  # Seconds per year
