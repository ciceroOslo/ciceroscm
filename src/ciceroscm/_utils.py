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
