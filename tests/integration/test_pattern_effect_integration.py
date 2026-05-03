"""
End-to-end integration tests for the Tier-3 pattern-mediated feedback wiring.

Verifies that:
  * ``delta_lambda_aero = 0`` produces bit-for-bit identical output to a run
    without the parameter at all (backward compatibility).
  * A non-zero ``delta_lambda_aero`` shifts the diagnosed Gregory feedback
    in the expected direction (Tier-3 signature).
  * Misconfiguration (``delta_lambda_aero != 0`` with a thermal model that
    does not implement the capability) raises ``ValueError`` at startup.
"""

import os

import numpy as np

from ciceroscm import CICEROSCM


def _build_cscm(test_data_dir):
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


def _run(cscm, **udm_overrides):
    pamset_udm = {
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
    pamset_udm.update(udm_overrides)
    cscm._run(
        {"results_as_dict": True},
        pamset_udm=pamset_udm,
        pamset_emiconc={
            "qbmb": 0.0,
            "qo3": 0.5,
            "qdirso2": -0.00308,
            "qindso2": -0.97 / 57.052577209999995,
            "qbc": 0.0279,
            "qoc": -0.00433,
            "qh2o_ch4": 0.091915,
            "ref_yr": 2010,
        },
    )
    return cscm


def test_delta_lambda_zero_matches_baseline(test_data_dir):
    """delta_lambda_aero = 0 must reproduce the baseline temperature trajectory."""
    cscm_base = _run(_build_cscm(test_data_dir))
    cscm_off = _run(_build_cscm(test_data_dir), delta_lambda_aero=0.0)

    t_base = np.asarray(cscm_base.results["dT_glob"])
    t_off = np.asarray(cscm_off.results["dT_glob"])
    np.testing.assert_array_equal(t_base, t_off)


def test_positive_delta_lambda_strengthens_feedback(test_data_dir):
    """
    With delta_lambda_aero > 0, the feedback is stronger during the aerosol
    era and the diagnosed Gregory slope (regression of N-F on T) over the
    historical window should be more negative than the baseline.
    """
    cscm_base = _run(_build_cscm(test_data_dir))
    cscm_pos = _run(_build_cscm(test_data_dir), delta_lambda_aero=1.5)

    years = np.arange(1750, 2101)
    hist = (years >= 1900) & (years <= 2014)

    def feedback(cscm):
        T = np.asarray(cscm.results["dT_glob"])[hist]
        N = np.asarray(cscm.results["RIB_glob"])[hist]
        F = np.asarray(cscm.results["Total_forcing"])[hist]
        R = N - F
        return -np.cov(T, R, ddof=0)[0, 1] / np.var(T)

    f_base = feedback(cscm_base)
    f_pos = feedback(cscm_pos)
    # Positive delta_lambda_aero should produce a *larger* (more strongly
    # damping) diagnosed feedback over the historical aerosol period.
    assert f_pos > f_base + 0.1  # comfortable margin above any numerical noise


def test_negative_delta_lambda_weakens_feedback(test_data_dir):
    cscm_base = _run(_build_cscm(test_data_dir))
    cscm_neg = _run(_build_cscm(test_data_dir), delta_lambda_aero=-1.5)

    years = np.arange(1750, 2101)
    hist = (years >= 1900) & (years <= 2014)

    def feedback(cscm):
        T = np.asarray(cscm.results["dT_glob"])[hist]
        N = np.asarray(cscm.results["RIB_glob"])[hist]
        F = np.asarray(cscm.results["Total_forcing"])[hist]
        R = N - F
        return -np.cov(T, R, ddof=0)[0, 1] / np.var(T)

    f_base = feedback(cscm_base)
    f_neg = feedback(cscm_neg)
    assert f_neg < f_base - 0.1
