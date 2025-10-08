"""
Functionality for carbon decay functions for both land (biotic/ rb) and ocean mixed layer (rs)
"""

import numpy as np


def _extract_indexed_parameters(pamset, prefix):
    """
    Extract indexed parameters from pamset.

    Parameters
    ----------
    pamset : dict
        Parameter set to extract from
    prefix : str
        Parameter prefix (e.g., "rb_coef", "rb_tim")

    Returns
    -------
    tuple
        (values_list, keys_to_remove)
    """
    values = []
    keys_to_remove = []

    for key in pamset.keys():
        if key.startswith(prefix) and key[len(prefix) :].isdigit():
            idx = int(key[len(prefix) :])
            # Ensure list is long enough
            while len(values) <= idx:
                values.append(None)
            values[idx] = pamset[key]
            keys_to_remove.append(key)

    # Remove None values
    values = [v for v in values if v is not None]
    return values, keys_to_remove


def _process_rb_parameters(pamset):
    """Process rb_function flat parameters."""
    rb_coefs, rb_coef_keys = _extract_indexed_parameters(pamset, "rb_coef")
    rb_tims, rb_tim_keys = _extract_indexed_parameters(pamset, "rb_tim")

    if not (rb_coefs or rb_tims):
        return pamset, []

    if len(rb_coefs) != len(rb_tims):
        raise ValueError(
            f"Number of rb_coef parameters ({len(rb_coefs)}) must match "
            f"number of rb_tim parameters ({len(rb_tims)})"
        )

    if rb_coefs:
        pamset["rb_function"] = {"coeffs": rb_coefs, "timescales": rb_tims}

    return pamset, rb_coef_keys + rb_tim_keys


def _process_rs_parameters(pamset):
    """Process rs_function flat parameters."""
    rs_coefs, rs_coef_keys = _extract_indexed_parameters(pamset, "rs_coef")
    rs_tims, rs_tim_keys = _extract_indexed_parameters(pamset, "rs_tim")

    if not (rs_coefs or rs_tims):
        return pamset, []

    # For rs_function, coeffs should have one more element than timescales
    if len(rs_coefs) != len(rs_tims) + 1:
        raise ValueError(
            f"For rs_function, number of rs_coef parameters ({len(rs_coefs)}) "
            f"must be one more than number of rs_tim parameters ({len(rs_tims)})"
        )

    if rs_coefs:
        pamset["rs_function"] = {"coeffs": rs_coefs, "timescales": rs_tims}

    return pamset, rs_coef_keys + rs_tim_keys


def _process_flat_carbon_parameters(pamset):
    """
    Process flat carbon cycle parameters and convert to dictionary format.

    This function allows users to specify carbon cycle function parameters
    as individual floats (e.g., rb_coef0, rb_coef1, rb_tim0, rb_tim1)
    instead of requiring dictionary structures.

    Parameters
    ----------
    pamset : dict
        Parameter set that may contain flat carbon cycle parameters

    Returns
    -------
    dict
        Updated pamset with flat parameters converted to dictionary structures

    Examples
    --------
    Input: {"rb_coef0": 0.5, "rb_coef1": 0.25, "rb_tim0": 2.5, "rb_tim1": 10}
    Output: {"rb_function": {"coeffs": [0.5, 0.25], "timescales": [2.5, 10]}}
    """
    # Early return if no flat carbon parameters are present
    has_flat_params = any(
        key.startswith(("rb_coef", "rb_tim", "rs_coef", "rs_tim"))
        for key in pamset.keys()
    )
    if not has_flat_params:
        return pamset

    pamset = pamset.copy()

    # Process rb_function parameters
    pamset, rb_keys_to_remove = _process_rb_parameters(pamset)

    # Process rs_function parameters
    pamset, rs_keys_to_remove = _process_rs_parameters(pamset)

    # Remove all flat parameter keys
    for key in rb_keys_to_remove + rs_keys_to_remove:
        pamset.pop(key)

    return pamset


def rb_function(it, idtm=24):
    """
    Calculate biotic decay function

    Calculate biotic decay function
    time is the year index*idtm + i, i.e. the month number

    Parameters
    ----------
    it : int
      is the time index, there are idtm time points per yer
    idtm : int
        Number of time points per year, default is 24
    Returns
    -------
    float
        The biotic decay function value for this time
    """
    time = it / idtm

    biotic_decay = (
        0.70211 * np.exp(-0.35 * time)
        + 13.4141e-3 * np.exp(-time / 20.0)
        - 0.71846 * np.exp(-55 * time / 120.0)
        + 2.9323e-3 * np.exp(-time / 100.0)
    )
    return biotic_decay


def rb_function2(it, rb_coef, rb_tim, idtm=24):
    """
    Calculate biotic decay function
    Calculate biotic decay function
    time is the year index*idtm + i, i.e. the month number

    Parameters
    ----------
    it : int
      is the time index, there are idtm time points per yer
    rb_coef : list
        Fractional attribution of pools (normalised to sum 1 in code)
        Example values could be [0.5,0.25,0.25]
    rb_tim : list
        Timescales in years for each pool
        example value [2.5,10,60]
    idtm : int
        Number of time points per year, default is 24


    Returns
    -------
    float
        The biotic decay function value for this time
        This function will have the following properties:
        It will integrate to idtm at infinity
        It will have value 0 at time 0
    """
    time = it / idtm
    if hasattr(time, "__len__"):
        nt = len(time)
    else:
        nt = 1
    ncoef = np.array(rb_coef) / np.sum(rb_coef)
    tcoef = ncoef / np.array(rb_tim) ** 2
    dmat = np.zeros((nt, len(rb_coef)))
    for i in range(len(rb_coef)):
        dmat[:, i] = (time * tcoef[i]).T * np.exp(-time / rb_tim[i])
    biotic_decay = np.sum(dmat, axis=1)
    if nt == 1:
        biotic_decay = biotic_decay[0]
    return biotic_decay


def rs_function(it, idtm=24):
    """
    Calculate pulse response function for mixed layer

    Calculate pulse response function for mixed layer
    time is the year index*idtm + i, i.e. the month number

    Parameters
    ----------
    it : int
      is the time index, there are idtm time points per yer
    idtm : int
        Number of time points per year, default is 24

    Returns
    -------
    float
         The pulse_response function for this time
    """
    time = it / idtm
    if time < 2.0:
        pulse_response = (
            0.12935
            + 0.21898 * np.exp(-time / 0.034569)
            + 0.17003 * np.exp(-time / 0.26936)
            + 0.24071 * np.exp(-time / 0.96083)
            + 0.24093 * np.exp(-time / 4.9792)
        )
    else:
        pulse_response = (
            0.022936
            + 0.24278 * np.exp(-time / 1.2679)
            + 0.13963 * np.exp(-time / 5.2528)
            + 0.089318 * np.exp(-time / 18.601)
            + 0.03782 * np.exp(-time / 68.736)
            + 0.035549 * np.exp(-time / 232.3)
        )
    return pulse_response


rs_function_array = np.vectorize(rs_function)


def rs_function2(it, rs_coef, rs_tim, idtm=24):
    """
    Calculate pulse response function for mixed layer

    Calculate pulse response function for mixed layer
    time is the year index*idtm + i, i.e. the month number

    Parameters
    ----------
    it : int
      is the time index, there are idtm time points per yer
    rs_coef : np.ndarray
        All the coefficients of the function, the first term
        is assumed to be a constant term, hence this array
        should have one more term, than the rs_tim array.
        In practice these will be normalised so they sum
        to one. Example values are [0.1, 0.6,0.15,0.15]
    rs_tim : np.ndarray
        All the timescales of the function. Should be one
        shorter than rs_coeff to account for a constant term.
        Example values are [.8,7,80]
    idtm : int
        Number of time points per year, default is 24

    Returns
    -------
    float
         The pulse_response function for this time
         This function has the following properites
         It should be 1 at time zero and asymptote
         towards a small non-negative value in the infinte future
    """
    time = it / idtm
    ncoef = np.array(rs_coef) / np.sum(rs_coef)
    if hasattr(time, "__len__"):
        nt = len(time)
    else:
        nt = 1
    dmat = np.zeros((len(rs_coef), nt))
    for i in range(len(rs_coef)):
        if i == 0:  # pylint: disable=compare-to-zero
            dmat[i, :] = ncoef[i]
        else:
            dmat[i, :] = ncoef[i] * np.exp(-time / rs_tim[i - 1])
    pulse_response = np.sum(dmat, axis=0)
    if nt == 1:
        pulse_response = pulse_response[0]
    return pulse_response
