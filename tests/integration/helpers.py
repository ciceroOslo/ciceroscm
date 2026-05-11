"""
Shared helpers for integration tests.

``build_cscm`` and ``run_cscm`` set up a standard SSP2-4.5 emissions run
with the default UDM and emiconc parameter sets used across the integration
suite. Tests that need a non-standard configuration should construct
CICEROSCM and call _run directly.
"""

import os

from ciceroscm import CICEROSCM

_DEFAULT_PAMSET_UDM = {
    "rlamdo": 16.0,
    "akapa": 0.634,
    "cpi": 0.4,
    "W": 4,
    "beto": 3.5,
    "lambda": 0.54,
    "mixed": 60.0,
    "foan": 0.61,
    "foas": 0.81,
    "ebbeta": 0.0,
    "fnso": 0.7531,
    "lm": 40,
    "ldtime": 12,
}

_DEFAULT_PAMSET_EMICONC = {
    "qbmb": 0.0,
    "qo3": 0.5,
    "qdirso2": -0.00308,
    "qindso2": -0.97 / 57.052577209999995,
    "qbc": 0.0279,
    "qoc": -0.00433,
    "qh2o_ch4": 0.091915,
    "ref_yr": 2010,
}


def build_cscm(test_data_dir):
    """Return an uninitialised CICEROSCM configured for the SSP2-4.5 test dataset."""
    return CICEROSCM(
        {
            "gaspam_file": os.path.join(test_data_dir, "gases_v1RCMIP.txt"),
            "nyend": 2100,
            "nystart": 1750,
            "emstart": 1850,
            "concentrations_file": os.path.join(test_data_dir, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(test_data_dir, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(test_data_dir, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(test_data_dir, "natemis_n2o.txt"),
        }
    )


def run_cscm(cscm, **udm_overrides):
    """Run *cscm* with the default parameter sets, applying any udm_overrides."""
    pamset_udm = {**_DEFAULT_PAMSET_UDM, **udm_overrides}
    cscm._run(
        {"results_as_dict": True},
        pamset_udm=pamset_udm,
        pamset_emiconc=_DEFAULT_PAMSET_EMICONC,
    )
    return cscm
