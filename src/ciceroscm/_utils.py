"""
Class for common utility functions
"""

import logging

LOGGER = logging.getLogger(__name__)


def check_numeric_pamset(required, pamset):
    """
    Check numeric pamset conforms

    Check that pamset has required inputs
    If not insert defaults from required dict

    Parameters
    ----------
    required : dict
            Dictionary containing the required keys
            with default values
    pamset : dict
          Dictionary to be checked and filled in from
          required where necessary

    Returns
    -------
    dict
         pamset augmented with necessary default values
         from required
    """
    for pam, value in required.items():
        if pam not in pamset:
            LOGGER.warning(  # pylint: disable=logging-fstring-interpolation
                f"Parameter {pam} not in pamset. Using default value {value}",
            )
            pamset[pam] = value
        elif not isinstance(pamset[pam], int) and not isinstance(pamset[pam], float):
            LOGGER.warning(  # pylint: disable=logging-fstring-interpolation
                f"Parameter {pam} must be a number. Using default value {value}",
            )
            pamset[pam] = value
    return pamset


def cut_non_required(required, pamset, cut_warnings=False):
    """
    Cut elements from pamset that are not required

    Take out elements of pamset that are not used

    Parameters
    ----------
    required : dict
            Dictionary containing the required keys
            with default values
    pamset : dict
          Dictionary to be checked and cut from
    cut_warnings : bool
                Boolean stating whether or not warnings
                should be printed if parameters are cut
    Returns
    -------
    dict
         pamset where unnecessary elements have been cut
    """
    new_pamset = {}
    for pam in pamset:
        if pam in required:
            new_pamset[pam] = pamset[pam]
        elif cut_warnings:
            LOGGER.warning(  # pylint: disable=logging-fstring-interpolation
                f"Parameter {pam} is not used. Please check if you have a typo",
            )
    return new_pamset


def cut_and_check_pamset(
    required, pamset, used={}, cut_warnings=False
):  # pylint: disable=dangerous-default-value
    """
    Combine cut and check of parameters

    Calling other methods to first cut unnecessary
    elements and afterwords put in missing elements

    Parameters
    ----------
    required : dict
            Dictionary containing the required keys
            with default values
    pamset : dict
          Dictionary to be checked and cut from
    used : dict
        Dictionary or list containing keys used, so
        that unused parameters can be cut.
    cut_warnings : bool
                Boolean stating whether or not warnings
                should be printed if parameters are cut
    Returns
    -------
    dict
         pamset where unnecessary elements have been cut
         and missing necessary elements have been added
         or augmented if not of numeric type.
    """
    if used:
        used.update(required)
    else:
        used = required
    pamset = cut_non_required(used, pamset, cut_warnings)
    return check_numeric_pamset(required, pamset)
