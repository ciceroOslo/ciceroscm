"""
Init
"""
from . import _version
from .ciceroscm import CICEROSCM

__version__ = _version.get_versions()["version"]
