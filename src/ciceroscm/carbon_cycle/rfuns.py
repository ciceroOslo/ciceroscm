"""
Functionality for carbon decay functions for both land (biotic/ rb) and ocean mixed layer (rs)
"""

import numpy as np


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
    rb_coefs = []
    rb_tims = []
    rb_keys_to_remove = []

    # Find all rb_coef* and rb_tim* parameters
    for key in pamset.keys():
        if key.startswith("rb_coef") and key[7:].isdigit():
            idx = int(key[7:])
            # Ensure list is long enough
            while len(rb_coefs) <= idx:
                rb_coefs.append(None)
            rb_coefs[idx] = pamset[key]
            rb_keys_to_remove.append(key)
        elif key.startswith("rb_tim") and key[6:].isdigit():
            idx = int(key[6:])
            # Ensure list is long enough
            while len(rb_tims) <= idx:
                rb_tims.append(None)
            rb_tims[idx] = pamset[key]
            rb_keys_to_remove.append(key)

    # Convert to rb_function dictionary if parameters found
    if rb_coefs or rb_tims:
        # Remove None values and check for consistency
        rb_coefs = [c for c in rb_coefs if c is not None]
        rb_tims = [t for t in rb_tims if t is not None]

        if len(rb_coefs) != len(rb_tims):
            raise ValueError(
                f"Number of rb_coef parameters ({len(rb_coefs)}) must match "
                f"number of rb_tim parameters ({len(rb_tims)})"
            )

        if len(rb_coefs) > 0:
            pamset["rb_function"] = {"coeffs": rb_coefs, "timescales": rb_tims}

        # Remove individual parameters
        for key in rb_keys_to_remove:
            pamset.pop(key)

    # Process rs_function parameters
    rs_coefs = []
    rs_tims = []
    rs_keys_to_remove = []

    # Find all rs_coef* and rs_tim* parameters
    for key in pamset.keys():
        if key.startswith("rs_coef") and key[7:].isdigit():
            idx = int(key[7:])
            # Ensure list is long enough
            while len(rs_coefs) <= idx:
                rs_coefs.append(None)
            rs_coefs[idx] = pamset[key]
            rs_keys_to_remove.append(key)
        elif key.startswith("rs_tim") and key[6:].isdigit():
            idx = int(key[6:])
            # Ensure list is long enough
            while len(rs_tims) <= idx:
                rs_tims.append(None)
            rs_tims[idx] = pamset[key]
            rs_keys_to_remove.append(key)

    # Convert to rs_function dictionary if parameters found
    if rs_coefs or rs_tims:
        # Remove None values and check for consistency
        rs_coefs = [c for c in rs_coefs if c is not None]
        rs_tims = [t for t in rs_tims if t is not None]

        # For rs_function, coeffs should have one more element than timescales
        if len(rs_coefs) != len(rs_tims) + 1:
            raise ValueError(
                f"For rs_function, number of rs_coef parameters ({len(rs_coefs)}) "
                f"must be one more than number of rs_tim parameters ({len(rs_tims)})"
            )

        if len(rs_coefs) > 0:
            pamset["rs_function"] = {"coeffs": rs_coefs, "timescales": rs_tims}

        # Remove individual parameters
        for key in rs_keys_to_remove:
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
