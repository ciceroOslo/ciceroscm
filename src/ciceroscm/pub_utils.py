"""
Module for public utilities

Presently mainly methods to make decay functions from arrays
"""

import logging

import numpy as np

LOGGER = logging.getLogger(__name__)


def _check_array_consistency(coeffs, timescales):
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
    if len(coeffs) != len(timescales):
        LOGGER.error("Coefficient array length must be equal to exponent array length")
        raise ValueError(
            "Coefficient array length must be equal to exponent array length"
        )
    if np.any(timescales < 0):
        LOGGER.error("All timescales must be positive numbers")
        raise ValueError("All timescales must be positive numbers")
    return coeffs, timescales


def make_rs_function_from_arrays(rs_coeff, rs_time):
    """
    Make mixed layer pulse response function from arrays of coefficients
    and timescales

    Parameters
    ----------
    rs_coeff : np.ndarray
        array of coefficients to the function
        These should all be postive and sum to less than one
        The constant term will be 1 minus this sum
    rs_time : np.ndarray
        array of decay timescales for each term
        All of these must be postive

    Returns
    -------
    function
       Mixed layer carbon pool decay function based on the coefficient and
       timescale vectors

    Raises
    ------
    ValueError
       If the coefficients are negative, or sum to more than one
    """
    rs_coeff, rs_time = _check_array_consistency(rs_coeff, rs_time)
    coeff0 = 1 - np.sum(rs_coeff)
    if np.any(rs_coeff <= 0):
        LOGGER.error("All coefficients must be positive")
        raise ValueError("All coefficients must be positive")
    if coeff0 < 0:
        raise ValueError(
            "Sum of carbon decay coefficients must be less than or equal to 1"
        )

    return lambda it, idtm: coeff0 + np.dot(rs_coeff, np.exp(-it / rs_time / idtm))


def make_rb_function_from_arrays(rb_coeff, rb_time):
    """
    Make carbon pool decay function from arrays of coefficients
    and timescales

    Parameters
    ----------
    rb_coeff : np.ndarray
        array of coefficients to the function
        All of these must be positive
    rb_time : np.ndarray
        array of decay timescales for each term
        All of these must be postive

    Returns
    -------
    function
       Biotic decay function based on the coefficient and
       timescales
    """
    rb_coeff, rb_time = _check_array_consistency(rb_coeff, rb_time)
    return lambda it, idtm: np.dot(rb_coeff, np.exp(-it / rb_time / idtm))


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
    print(numpart)
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
