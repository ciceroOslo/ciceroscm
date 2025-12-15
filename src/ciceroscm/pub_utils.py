"""
Module for public utilities

Presently mainly methods to make decay functions from arrays
"""

import logging

import numpy as np

LOGGER = logging.getLogger(__name__)

REG_MAPPING_DEFAULT = {
    "ASIA": "EAS",
    "LAM": "SAM",
    "MAF": "NAF",
    "OECD": "NAM",
    "REF": "RBU",
}


def get_first_key(dct, value):
    """
    Return the first key in a dictionary that matches the given value.

    Parameters
    ----------
    dct : dict
        The dictionary to search.
    value : any
        The value to search for in the dictionary.

    Returns
    -------
    key : any
        The first key in the dictionary whose value matches the specified value.

    Raises
    ------
    IndexError
        If the value is not found in the dictionary.
    """
    return [key for key in dct if (dct[key] == value)][0]


def _check_array_consistency(coeffs, timescales, for_rs=False):
    """
    Check array consistency for sets of coefficients and
    decay timescales

    Parameters
    ----------
    coeffs : array
        array, preferably an np.array
    timescales : array
        array, preferably an np.array

    Returns
    -------
    list
        list of two np.ndarrays that are coeffs and timescales
        converted to numpy arrays in case they weren't before
        and checked to be of equal length and with the timescales
        guaranteed to be positive in all errors

    Raises
    ------
    TypeError
        If coeffs or timescales are not arrays
    ValueError
        If any of the timescales are negative or
        if length of timescales and coeffs don't match
    """
    try:
        coeffs = np.array(coeffs)
        timescales = np.array(timescales)
    except TypeError as terr:
        LOGGER.error("Both coefficients and timescales must be 1D arrays")
        raise terr
    if for_rs and len(coeffs) != (len(timescales) + 1):
        LOGGER.error(
            "Coefficient array must have one more element than timescale array"
        )
        raise ValueError(
            "Coefficient array length must be equal to exponent array length + 1"
        )
    if not for_rs and len(coeffs) != len(timescales):
        LOGGER.error("Coefficient array must be the same length as the timescale array")
        raise ValueError(
            "Coefficient array length must be equal to exponent array length"
        )
    if np.any(timescales < 0):
        LOGGER.error("All timescales must be positive numbers")
        raise ValueError("All timescales must be positive numbers")
    if np.any(timescales < 0):
        LOGGER.error("All timescales must be positive numbers")
        raise ValueError("All timescales must be positive numbers")
    if np.any(coeffs < 0):
        LOGGER.error("All coefficients must be positive numbers")
        raise ValueError("All coefficients must be positive numbers")
    return coeffs, timescales


def find_num_cl_in_hcfc(tracer):
    """
    Find number of Cl-atoms in a CFC or HCFC tracer

    Assuming for now, that no bromines are involved in these
    and Halon names are used for bromine containing CFCs

    Parameters
    ----------
    tracer : str
             Name of tracer, should start with HCFC or CFC
    Returns
    -------
    int
       Number of chlorine atoms in the compound
    """
    numpart = tracer.split("-")[-1]
    if not numpart[-1].isdigit():
        numpart = int(numpart[:-1])
    else:
        numpart = int(numpart)
    double_bonds = numpart // 1000
    bonds_total = ((numpart % 1000) // 100) * 2 + 4
    hydrogen = (numpart % 100) // 10 - 1
    fluorine = numpart % 10
    chlorine = bonds_total - 2 * double_bonds - hydrogen - fluorine
    return chlorine


def make_cl_and_br_dictionaries(gases_index):
    """
    Make dictionaries with compound names as keys
    and number of bromine or chlorine as values
    from the compound names in a list
    (from gaspamfile headers)

    Parameters
    ----------
    gases_index : list
                 List with compound names as strings
                 This should be the index of the
                 DataFrame generated from a gaspam_file
    Returns
    -------
    tuple
         The chlorine dictionary followed by the
         bromine dictionary filled in
    """
    chlor_dict = {}
    brom_dict = {}
    for tracer in gases_index:
        if "CFC-" in tracer:
            chlor_dict[tracer] = find_num_cl_in_hcfc(tracer)
        elif tracer[:2] == "H-":
            brom_dict[tracer] = int(tracer[5])
            if int(tracer[4]) > 0:
                chlor_dict[tracer] = int(tracer[4])
        elif "Cl" in tracer:
            if len(tracer.split("Cl")[-1]) > 0 and tracer.split("Cl")[-1][0].isdigit():
                chlor_dict[tracer] = int(tracer.split("Cl")[-1][0])
            else:
                chlor_dict[tracer] = 1
        elif "Br" in tracer:
            if len(tracer.split("Br")[-1]) > 0 and tracer.split("Br")[-1][0].isdigit():
                brom_dict[tracer] = int(tracer.split("Br")[-1][0])
            else:
                brom_dict[tracer] = 1

    return chlor_dict, brom_dict


def make_regional_aerosol_gaspamdata(
    gases_data,
    reg_aerosol_rf_data,
    reg_mapping=None,
    unit_conversion=1e9,  # Convert W/m² per Tg to W/m² per kg (1e9 kg = 1 Tg)
):
    """
    Make regional aerosol gaspamdata from existing gaspam and some regional data

    Creates separate gas parameter entries for each aerosol-region combination
    (e.g., SO2_ASIA, BC_OECD) using region-specific forcing efficacies.

    Parameters
    ----------
    gases_data : pd.DataFrame
        Original gaspamdata without regional aerosols
    reg_aerosol_rf_data : pd.DataFrame
        Dataframe with regional forcing for aerosols (rows=regions, cols=aerosols)
    reg_mapping : dict, optional
        Mapping between emission regions and forcing regions
        Default maps ASIA→EAS, LAM→SAM, MAF→NAF, OECD→NAM, REF→RBU
    unit_conversion : float, optional
        Unit conversion factor from forcing data units to model units (default: 1e9)

    Returns
    -------
    pd.DataFrame
        New gaspam dataframe with extra entries for the added regional aerosols
    """
    if reg_mapping is None:
        reg_mapping = REG_MAPPING_DEFAULT
    gases_new = gases_data.copy()
    for column in reg_aerosol_rf_data.columns:
        row_start = gases_data.loc[column].copy()
        for key, expression in reg_mapping.items():
            row_start["ALPHA"] = (
                reg_aerosol_rf_data[column][expression] * unit_conversion
            )
            gases_new.loc[f"{column}_{key}"] = row_start
    return gases_new
