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

# Pattern-mediated feedback: tracers grouped as "aerosol" for the
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

# Fractional release factors used for equivalent effective stratospheric
# chlorine (EESC) calculations, per halocarbon species.
EESC_FRF = {
    "CFC-11": 0.47,
    "CFC-12": 0.24,
    "CFC-13": 0.06,
    "CFC-112": 0.3,
    "CFC-112a": 0.53,
    "CFC-113": 0.30,
    "CFC-113a": 0.29,
    "CFC-114": 0.13,
    "CFC-114a": 0.32,
    "CFC-115": 0.07,
    "CCl4": 0.56,
    "CH3CCl3": 0.61,
    "HCFC-133a": 0.4,
    "HCFC-22": 0.15,
    "HCFC-141b": 0.34,
    "HCFC-142b": 0.17,
    "HCFC-123": 0.66,
    "HCFC-124": 0.32,
    "HCFC-31": 0.67,
    "H-1211": 0.65,
    "H-1202": 0.67,
    "H-1301": 0.32,
    "H-2402": 0.66,
    "CH3Br": 0.6,
    "CH3Cl": 0.44,
}

# Age-of-air weighting used to smear yearly EESC source terms over the
# preceding 20 years (age-spectrum convolution kernel).
EESC_WEIGHTS = [
    2.51540874e-01,
    2.85588337e-01,
    1.83647911e-01,
    1.09745490e-01,
    6.53733900e-02,
    3.93816364e-02,
    2.40508141e-02,
    1.48813664e-02,
    9.31632867e-03,
    5.89277991e-03,
    3.76106432e-03,
    2.41954343e-03,
    1.56739065e-03,
    1.02163109e-03,
    6.69558679e-04,
    4.40970373e-04,
    2.91703198e-04,
    1.93730846e-04,
    1.29128440e-04,
    8.63517467e-05,
]
