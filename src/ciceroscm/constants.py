"""
Physical and mathematical constants used throughout CICERO-SCM

This module contains fundamental physical constants, conversion factors,
and commonly used values across the climate model components.

Constants are organized by category and include detailed documentation
about their usage and derivation where applicable.
"""

# Time conversion constants
SEC_DAY = 86400  # Seconds per day
DAY_YEAR = 365.0  # Days per year (simplified, no leap years)

# Water properties (seawater)
WATER_DENSITY = 1000  # kg/m^3 - Density of seawater
WATER_HEAT_CAPACITY = 4181  # J/kg/K - Specific heat capacity of seawater

# Ocean model specific constants (UDM)
SEAWATER_DENSITY_UDM = 1030.0  # kg/m^3 - Seawater density used in UDM OHC calculations
SEAWATER_DENSITY_OHC = (
    1.03  # g/cm^3 - Density for UDM heat capacity calculations (legacy value)
)
SEAWATER_HEAT_CAPACITY_UDM = 0.955  # Dimensionless heat capacity factor for UDM
UDM_CONVERSION_FACTOR = 0.485  # Dimensionless conversion factor for UDM heat capacity
UDM_OHC_CONSTANT = 3.997e-19  # Conversion constant for ocean heat content calculations (units: see UDM docs)
HEMISPHERE_OCEAN_AREA = 2.55e14  # m^2 - Ocean area per hemisphere (used in UDM)
OCEAN_LAYER_THICKNESS = 100.0  # m - Standard ocean layer thickness in UDM
DEPTH_700M_LAYER_INDEX = 7  # Layer index constant for 700m depth calculations in UDM

# Earth surface areas
OCEAN_AREA = 3.61e14  # m^2 - Global ocean surface area (used in Two-Layer model)

# Derived constants
YEAR_IN_SECONDS = DAY_YEAR * SEC_DAY  # Seconds per year
