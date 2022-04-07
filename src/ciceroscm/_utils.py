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


def cut_non_required(required, pamset):
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

    Returns
    -------
    dict
         pamset where unnecessary elements have been cut
    """
    new_pamset = {}
    for pam in required:
        if pam in pamset:
            new_pamset[pam] = pamset[pam]
    return new_pamset


def cut_and_check_pamset(required, pamset):
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

    Returns
    -------
    dict
         pamset where unnecessary elements have been cut
         and missing necessary elements have been added
         or augmented if not of numeric type.
    """
    pamset = cut_non_required(required, pamset)
    return check_numeric_pamset(required, pamset)
