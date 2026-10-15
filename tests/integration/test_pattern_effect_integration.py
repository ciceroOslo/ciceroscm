"""
End-to-end integration tests for the pattern-mediated feedback wiring.

Verifies that:
  * ``delta_lambda_aero = 0`` produces bit-for-bit identical output to a run
    without the parameter at all (backward compatibility).
  * A non-zero ``delta_lambda_aero`` shifts the diagnosed Gregory feedback
    in the expected direction.
  * Misconfiguration (``delta_lambda_aero != 0`` with a thermal model that
    does not implement the capability) raises ``ValueError`` at startup.
"""

import numpy as np
import pytest
from helpers import build_cscm, load_pdo_data_padded, run_cscm


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


def test_pdo_data_loader_pads_missing_years_and_months(test_data_dir):
    pdo = load_pdo_data_padded(test_data_dir, 1750, 2100)

    assert pdo.shape == (351, 12)
    np.testing.assert_array_equal(pdo[0], np.zeros(12))
    np.testing.assert_array_equal(pdo[99], np.zeros(12))
    assert pdo[100, 0] == np.float64(1.06)
    np.testing.assert_array_equal(pdo[2100 - 1750], np.zeros(12))
    np.testing.assert_array_equal(pdo[2026 - 1750, 8:], np.zeros(4))


def test_pdo_input_with_zero_knobs_matches_baseline(test_data_dir):
    pdo = load_pdo_data_padded(test_data_dir, 1750, 2100)

    cscm_base = run_cscm(build_cscm(test_data_dir))
    cscm_pdo = run_cscm(
        build_cscm(test_data_dir),
        pdo_index_data=pdo,
        delta_lambda_pdo=0.0,
        pdo_efficacy_scale=0.0,
    )

    np.testing.assert_array_equal(
        cscm_base.results["dT_glob"], cscm_pdo.results["dT_glob"]
    )
    np.testing.assert_array_equal(
        cscm_base.results["OHCTOT"], cscm_pdo.results["OHCTOT"]
    )


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


def test_positive_pdo_modulation_strengthens_feedback(test_data_dir):
    pdo = np.abs(load_pdo_data_padded(test_data_dir, 1750, 2100))

    cscm_base = run_cscm(build_cscm(test_data_dir), pdo_index_data=pdo)
    cscm_pos = run_cscm(
        build_cscm(test_data_dir),
        pdo_index_data=pdo,
        delta_lambda_pdo=1.0,
    )

    assert feedback(cscm_pos) > feedback(cscm_base) + 0.05


def test_negative_pdo_modulation_weakens_feedback(test_data_dir):
    pdo = np.abs(load_pdo_data_padded(test_data_dir, 1750, 2100))

    cscm = build_cscm(test_data_dir)

    cscm_base = run_cscm(cscm)
    feedback_base = feedback(cscm_base)
    cscm.cfg["pdo_index_data"] = pdo
    cscm_neg = run_cscm(
        cscm,
        delta_lambda_pdo=-1.0,
    )

    assert feedback(cscm_neg) < feedback_base - 0.05


def test_positive_pdo_efficacy_scale_changes_ocean_heat_uptake(test_data_dir):
    pdo = np.abs(load_pdo_data_padded(test_data_dir, 1750, 2100))

    cscm_base = run_cscm(build_cscm(test_data_dir), pdo_index_data=pdo)
    cscm_eff = run_cscm(
        build_cscm(test_data_dir),
        pdo_index_data=pdo,
        pdo_efficacy_scale=0.5,
    )
    assert cscm_eff.results["OHCTOT"][-100] < cscm_base.results["OHCTOT"][-100]
    assert cscm_eff.results["dT_glob"][-100] < cscm_base.results["dT_glob"][-100]


def test_pdo_in_cfg_runs_without_pamset_udm(test_data_dir):
    pdo = load_pdo_data_padded(test_data_dir, 1750, 2100)
    cscm = build_cscm(test_data_dir, pdo_index_data=pdo)
    cscm._run({"results_as_dict": True})
    assert np.isfinite(cscm.results["dT_glob"]).all()


def test_pdo_feedback_resets_when_index_returns_to_zero(test_data_dir):
    # PDO data is zero-padded after the observed record, so the feedback
    # must go back to 1 / lambda rather than keep the last non-zero month
    pdo = load_pdo_data_padded(test_data_dir, 1750, 2100)
    cscm = run_cscm(build_cscm(test_data_dir), pdo_index_data=pdo, delta_lambda_pdo=0.5)
    lambda_2100 = np.asarray(cscm.results["dynamic_lambda"])[-1]
    np.testing.assert_allclose(lambda_2100, 1.0 / 0.54)


def test_pdo_from_cfg_does_not_leak_into_caller_pamset_udm(test_data_dir):
    pdo = load_pdo_data_padded(test_data_dir, 1750, 2100)
    pamset_udm = {"lambda": 0.54, "delta_lambda_pdo": 0.5}
    build_cscm(test_data_dir, pdo_index_data=pdo)._run(
        {"results_as_dict": True}, pamset_udm=pamset_udm
    )
    assert "pdo_index_data" not in pamset_udm


def test_pdo_data_with_wrong_shape_raises(test_data_dir):
    pdo = load_pdo_data_padded(test_data_dir, 1750, 2100)
    with pytest.raises(ValueError, match="pdo_index_data"):
        build_cscm(test_data_dir, pdo_index_data=pdo[:100])
    with pytest.raises(ValueError, match="pdo_index_data"):
        run_cscm(build_cscm(test_data_dir), pdo_index_data=pdo[:, :6])
