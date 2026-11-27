import numpy as np
import pytest
from scipy.linalg.lapack import get_lapack_funcs

from ciceroscm.thermal_model import _udm_kernel
from ciceroscm.thermal_model.upwelling_diffusion_model import UpwellingDiffusionModel

GTSV = get_lapack_funcs("gtsv", (np.empty(1),))


@pytest.mark.parametrize("dominant", [True, False])
def test_gtsv_matches_lapack_bitwise(dominant):
    rng = np.random.default_rng(0)
    for _ in range(200):
        n = 40
        a, c = rng.standard_normal(n), rng.standard_normal(n)
        b = rng.standard_normal(n) + (4.0 if dominant else 0.0)
        d = rng.standard_normal(n)
        expected = GTSV(a[:-1], b, c[1:], d)[3]
        x = d.copy()
        _udm_kernel._gtsv(a.copy(), b.copy(), c.copy(), x)
        np.testing.assert_array_equal(x, expected)


@pytest.mark.parametrize("n", [1, 7, 12, 39, 40, 300])
def test_sum_matches_numpy_bitwise(n):
    rng = np.random.default_rng(n)
    for _ in range(200):
        x = rng.standard_normal(n) * 10.0 ** rng.integers(-8, 3, n)
        assert _udm_kernel._sum(x) == np.sum(x)


@pytest.mark.parametrize("efficacy", [1.0, 1.4])
def test_coeff_matches_model_bitwise(efficacy):
    udm = UpwellingDiffusionModel({"ocean_efficacy": efficacy, "ebbeta": 0.3})
    udm.setup_ebud2(0.4, 0.2)
    wcfac = udm.pamset["W"] / (86400 * 365.0) * (1 - 0.3 * 0.4 / 7.0)
    wcfac = wcfac * udm.pamset["dt"]
    abc = np.empty((3, udm.pamset["lm"]))
    gff = udm.get_gam_and_fro_factor_ns(True)
    _udm_kernel._coeff(abc, wcfac, gff, efficacy, udm._kernel[1][3])
    for i, k in enumerate("abc"):
        np.testing.assert_array_equal(abc[i], udm.varrying[f"{k}coeffn"])


def _forcing(nyears):
    rng = np.random.default_rng(1)
    forc = np.cumsum(rng.normal(0.05, 0.3, (nyears, 2)), axis=0)
    volc = np.where(rng.random((nyears, 1)) < 0.1, -rng.random((nyears, 12)), 0.0)
    return forc, volc, rng.random(nyears)


PAMSETS = [
    {},
    {"threstemp": 0, "ebbeta": 0.5, "ocean_efficacy": 1.3},
    {"delta_lambda_aero": 0.3, "lm": 30, "ldtime": 24},
]


@pytest.mark.parametrize("pamset", PAMSETS)
def test_run_years_matches_energy_budget(pamset):
    forc, volc, w_aero = _forcing(30)
    by_year = UpwellingDiffusionModel(dict(pamset))
    res = [
        by_year.energy_budget([*forc[i], volc[i], volc[i]], w_aero[i], i)
        for i in range(len(forc))
    ]
    udm = UpwellingDiffusionModel(dict(pamset))
    together = udm.run_years(forc, volc, volc, w_aero)
    for key, value in together.items():
        np.testing.assert_array_equal(value, [r[key] for r in res])
    np.testing.assert_array_equal(udm.tn, by_year.tn)


def test_pdo_feedback_changes_within_year():
    pdo = np.tile(np.linspace(-1, 1, 12), (5, 1))
    udm = UpwellingDiffusionModel(
        {"pdo_index_data": pdo, "delta_lambda_pdo": 0.2, "pdo_efficacy_scale": 0.1}
    )
    forc, volc, _ = _forcing(5)
    res = udm.run_years(forc, volc, volc)
    base = 1.0 / udm.pamset["lambda"]
    assert res["dynamic_lambda"][0] == pytest.approx(base)
    assert udm.pamset["rlamda"] == pytest.approx(base + 0.2)
    assert res["pdo_index"][-1] == 1.0


@pytest.mark.skipif(
    not hasattr(_udm_kernel.run_years, "py_func"), reason="numba not installed"
)
def test_compiled_matches_python():
    forc, volc, _ = _forcing(20)
    udm = UpwellingDiffusionModel({"ocean_efficacy": 1.2})
    rlamda = np.full((20, 12), udm.pamset["rlamda"])
    efficacy = np.full((20, 12), 1.2)
    out = [np.empty((20, len(_udm_kernel.OUT_KEYS))) for _ in range(2)]
    for run, o in zip((_udm_kernel.run_years, _udm_kernel.run_years.py_func), out):
        state = (np.zeros(40), np.zeros(40), np.zeros(2))
        run(state, forc, volc, volc, rlamda, efficacy, *udm._kernel, o)
    np.testing.assert_array_equal(*out)
