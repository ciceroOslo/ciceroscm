"""
Class for common utility functions
"""
import logging

LOGGER = logging.getLogger(__name__)


def check_numeric_pamset(required, pamset):
    """
    Check that pamset has required inputs
    If not insert defaults from required dict
    Returns updated pamset
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
    Returns updated pamset
    """
    new_pamset = {}
    for pam in required:
        if pam in pamset:
            new_pamset[pam] = pamset[pam]
    return new_pamset


def cut_and_check_pamset(required, pamset):
    """
    Combine cutting pamset to required subset and
    checking the resulting subset has required values
    Returns updated pamset
    """
    pamset = cut_non_required(required, pamset)
    return check_numeric_pamset(required, pamset)
