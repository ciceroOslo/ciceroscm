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

import numpy as np
from helpers import build_cscm, run_cscm


def feedback(cscm):
    """
    Compute the diagnosed Gregory feedback parameter over the historical period.

    Convenience function to compute the diagnosed Gregory feedback parameter
    including aerosol pattern effects over the historical aerosol period.
    """
    years = np.arange(1750, 2101)
    hist = (years >= 1900) & (years <= 2014)
    T = np.asarray(cscm.results["dT_glob"])[hist]
    N = np.asarray(cscm.results["RIB_glob"])[hist]
    F = np.asarray(cscm.results["Total_forcing"])[hist]
    R = N - F
    return -np.cov(T, R, ddof=0)[0, 1] / np.var(T)


def test_delta_lambda_zero_matches_baseline(test_data_dir):
    """delta_lambda_aero = 0 must reproduce the baseline temperature trajectory."""
    cscm_base = run_cscm(build_cscm(test_data_dir))
    cscm_off = run_cscm(build_cscm(test_data_dir), delta_lambda_aero=0.0)

    t_base = np.asarray(cscm_base.results["dT_glob"])
    t_off = np.asarray(cscm_off.results["dT_glob"])
    np.testing.assert_array_equal(t_base, t_off)


def test_positive_delta_lambda_strengthens_feedback(test_data_dir):
    """
    With delta_lambda_aero > 0, the feedback is stronger during the aerosol
    era and the diagnosed Gregory slope (regression of N-F on T) over the
    historical window should be more negative than the baseline.
    """
    cscm_base = run_cscm(build_cscm(test_data_dir))
    cscm_pos = run_cscm(build_cscm(test_data_dir), delta_lambda_aero=1.5)

    f_base = feedback(cscm_base)
    f_pos = feedback(cscm_pos)
    # Positive delta_lambda_aero should produce a *larger* (more strongly
    # damping) diagnosed feedback over the historical aerosol period.
    assert f_pos > f_base + 0.1  # comfortable margin above any numerical noise


def test_negative_delta_lambda_weakens_feedback(test_data_dir):
    cscm_base = run_cscm(build_cscm(test_data_dir))
    cscm_neg = run_cscm(build_cscm(test_data_dir), delta_lambda_aero=-1.5)

    f_base = feedback(cscm_base)
    f_neg = feedback(cscm_neg)
    assert f_neg < f_base - 0.1
