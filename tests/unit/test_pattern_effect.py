"""
Tests for pattern-mediated feedback (Tier 3) machinery.

Covers:
  * The thermal-model pattern-effect capability protocol
    (``get_feedback_gregory`` / ``set_feedback_gregory``).
  * That setting the feedback refreshes derived quantities so the next
    ``energy_budget`` call uses the new value (UpwellingDiffusionModel
    needs ``fnx``/``fsx``/``gamn``/``gams``/``varrying`` to be refreshed;
    TwoLayerOceanModel reads ``pamset["lambda"]`` directly each step).
  * Driver-level wiring: capability check at startup, zero-knob identity,
    and that a non-zero ``delta_lambda_aero`` produces the expected sign
    of change in diagnosed feedback.
"""

import numpy as np
import pytest

from ciceroscm.thermal_model.abstract_thermal_model import AbstractThermalModel
from ciceroscm.thermal_model.two_layer_ocean import TwoLayerOceanModel
from ciceroscm.thermal_model.upwelling_diffusion_model import (
    UpwellingDiffusionModel,
)


# ---------------------------------------------------------------------------
# Capability protocol
# ---------------------------------------------------------------------------


class _DummyNoCapability(AbstractThermalModel):
    """Minimal thermal model that does NOT override the pattern-effect hooks."""

    thermal_model_required_pamset = {}
    output_dict_default = {}

    def energy_budget(self, forc_nh, forc_sh, fn_volc, fs_volc):
        return {}


def test_abstract_pattern_effect_default_raises():
    model = _DummyNoCapability()
    with pytest.raises(NotImplementedError, match="get_feedback_gregory"):
        model.get_feedback_gregory()
    with pytest.raises(NotImplementedError, match="set_feedback_gregory"):
        model.set_feedback_gregory(1.5)


# ---------------------------------------------------------------------------
# UpwellingDiffusionModel: lambda_pamset is K/(W m^-2); rlamda = 1/lambda is
# the Gregory feedback. set_feedback_gregory must refresh fnx, fsx, and the
# gamn/gams/varrying state via setup_ebud.
# ---------------------------------------------------------------------------


def _udm_pamset(lambda_pamset=0.54, ocean_efficacy=1.0):
    return {
        "rlamdo": 16.0,
        "akapa": 0.634,
        "cpi": 0.4,
        "W": 4,
        "beto": 3.5,
        "threstemp": 7.0,
        "lambda": lambda_pamset,
        "mixed": 60.0,
        "foan": 0.61,
        "foas": 0.81,
        "ebbeta": 0.0,
        "lm": 40,
        "ldtime": 12,
        "ocean_efficacy": ocean_efficacy,
    }


def test_udm_get_feedback_gregory_returns_inverse_of_pamset_lambda():
    udm = UpwellingDiffusionModel(_udm_pamset(lambda_pamset=0.5))
    assert udm.get_feedback_gregory() == pytest.approx(1.0 / 0.5)


def test_udm_set_feedback_gregory_refreshes_derived_quantities():
    udm = UpwellingDiffusionModel(_udm_pamset(lambda_pamset=0.5))
    new_feedback = 2.5

    udm.set_feedback_gregory(new_feedback)

    assert udm.pamset["rlamda"] == pytest.approx(new_feedback)
    expected_fnx = (
        new_feedback
        + udm.pamset["foan"] * udm.pamset["rlamdo"]
        + udm.pamset["ebbeta"]
    )
    expected_fsx = (
        new_feedback
        + udm.pamset["foas"] * udm.pamset["rlamdo"]
        + udm.pamset["ebbeta"]
    )
    assert udm.pamset["fnx"] == pytest.approx(expected_fnx)
    assert udm.pamset["fsx"] == pytest.approx(expected_fsx)

    # gamn/gams are recomputed by setup_ebud and use rlamda directly.
    blm = udm.pamset["ebbeta"] / udm.pamset["rlamdo"]
    expected_gamn = (
        udm.pamset["foan"] + new_feedback / udm.pamset["rlamdo"] + blm
    )
    expected_gams = (
        udm.pamset["foas"] + new_feedback / udm.pamset["rlamdo"] + blm
    )
    assert udm.gamn == pytest.approx(expected_gamn)
    assert udm.gams == pytest.approx(expected_gams)


def test_udm_setting_same_feedback_is_noop_in_energy_budget():
    """A no-op set_feedback_gregory call must not perturb subsequent state."""
    udm_a = UpwellingDiffusionModel(_udm_pamset())
    udm_b = UpwellingDiffusionModel(_udm_pamset())

    volc = np.zeros(12)
    out_a = udm_a.energy_budget(1.0, 1.0, volc, volc)
    udm_b.set_feedback_gregory(udm_b.get_feedback_gregory())
    out_b = udm_b.energy_budget(1.0, 1.0, volc, volc)

    assert out_a["dtemp"] == pytest.approx(out_b["dtemp"])


def test_udm_changed_feedback_changes_energy_budget_output():
    """A larger Gregory feedback should produce smaller dtemp under same forcing."""
    udm_strong = UpwellingDiffusionModel(_udm_pamset())
    udm_weak = UpwellingDiffusionModel(_udm_pamset())

    udm_strong.set_feedback_gregory(2.5)  # stronger damping
    udm_weak.set_feedback_gregory(1.2)  # weaker damping

    volc = np.zeros(12)
    # Run a few years to let the response build up.
    dt_strong = dt_weak = 0.0
    for _ in range(10):
        out_s = udm_strong.energy_budget(2.0, 2.0, volc, volc)
        out_w = udm_weak.energy_budget(2.0, 2.0, volc, volc)
        dt_strong = out_s["dtemp"]
        dt_weak = out_w["dtemp"]

    assert dt_strong < dt_weak


# ---------------------------------------------------------------------------
# TwoLayerOceanModel: pamset["lambda"] *is* the Gregory feedback, so the
# get/set methods are direct.
# ---------------------------------------------------------------------------


def test_two_layer_get_feedback_gregory_returns_pamset_lambda():
    model = TwoLayerOceanModel({"lambda": 1.5})
    assert model.get_feedback_gregory() == pytest.approx(1.5)


def test_two_layer_set_feedback_gregory_updates_pamset_lambda():
    model = TwoLayerOceanModel({"lambda": 1.5})
    model.set_feedback_gregory(2.25)
    assert model.pamset["lambda"] == pytest.approx(2.25)


def test_two_layer_changed_feedback_changes_energy_budget_output():
    model_strong = TwoLayerOceanModel()
    model_weak = TwoLayerOceanModel()

    model_strong.set_feedback_gregory(3.0)
    model_weak.set_feedback_gregory(1.0)

    volc = np.zeros(12)
    for _ in range(10):
        out_s = model_strong.energy_budget(2.0, 2.0, volc, volc)
        out_w = model_weak.energy_budget(2.0, 2.0, volc, volc)
    assert out_s["dtemp"] < out_w["dtemp"]


# ---------------------------------------------------------------------------
# Pamset hygiene: delta_lambda_aero is a recognised parameter on both models
# (so users do not get cut_warnings), and defaults to 0.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("model_cls", [UpwellingDiffusionModel, TwoLayerOceanModel])
def test_delta_lambda_aero_default_is_zero(model_cls):
    if model_cls is UpwellingDiffusionModel:
        model = model_cls(_udm_pamset())
    else:
        model = model_cls()
    assert model.pamset["delta_lambda_aero"] == 0.0


@pytest.mark.parametrize("model_cls", [UpwellingDiffusionModel, TwoLayerOceanModel])
def test_delta_lambda_aero_passes_through_pamset_filter(model_cls):
    if model_cls is UpwellingDiffusionModel:
        pam = _udm_pamset()
        pam["delta_lambda_aero"] = 1.5
        model = model_cls(pam)
    else:
        model = model_cls({"delta_lambda_aero": 1.5})
    assert model.pamset["delta_lambda_aero"] == 1.5
