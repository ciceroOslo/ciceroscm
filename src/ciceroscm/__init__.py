"""
Init
"""

from . import _version
from .ciceroscm import CICEROSCM  # noqa: F401

__version__ = _version.get_versions()["version"]
