"""
Compiled year loop of the upwelling-diffusion energy budget

``run_years`` steps the ocean through one or more years. It performs the
same floating-point operations in the same order as the original numpy
implementation, including LAPACK's tridiagonal elimination and numpy's
pairwise summation, so results are bit-identical. It is compiled with
numba when available and otherwise runs as plain Python.
"""

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover

    def njit(*args, **kwargs):  # pylint: disable=unused-argument
        """Run undecorated when numba is not installed."""
        return lambda f: f


_JIT = njit(cache=True, error_model="numpy")

# Indices of the per-year results written by run_years
OUT_KEYS = (
    "dtempnh",
    "dtempsh",
    "dtempnh_air",
    "dtempsh_air",
    "dtempnh_sea",
    "dtempsh_sea",
    "RIBN",
    "RIBS",
    "anomalous_radiation",
    "OHC700",
    "OHCTOT",
)


@_JIT
def _sum(x):
    """Sum ``x`` in numpy's pairwise order (bit-identical to ``np.sum``)."""
    n = len(x)
    if n < 8:  # noqa: PLR2004
        res = 0.0
        for i in range(n):
            res += x[i]
        return res
    if n <= 128:  # noqa: PLR2004
        r = x[:8].copy()
        m = n - n % 8
        for i in range(8, m, 8):
            for j in range(8):
                r[j] += x[i + j]
        res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]))
        for i in range(m, n):
            res += x[i]
        return res
    n2 = n // 2
    n2 -= n2 % 8
    return _sum(x[:n2]) + _sum(x[n2:])


@_JIT
def _gtsv(a, b, c, x):
    """
    Solve the tridiagonal system in place as LAPACK ``dgtsv`` does

    Row i reads ``a[i-1] x[i-1] + b[i] x[i] + c[i+1] x[i+1] = d[i]``;
    ``x`` holds ``d`` on entry and the solution on exit, and ``a``, ``b``
    and ``c`` are overwritten.
    """
    n = len(x)
    dl, dd, du = a, b, c[1:]
    for i in range(n - 1):
        if abs(dd[i]) >= abs(dl[i]):
            fact = dl[i] / dd[i]
            dd[i + 1] = dd[i + 1] - fact * du[i]
            x[i + 1] = x[i + 1] - fact * x[i]
            dl[i] = 0.0
        else:  # partial pivoting: interchange rows i and i+1
            fact = dd[i] / dl[i]
            dd[i] = dl[i]
            temp = dd[i + 1]
            dd[i + 1] = du[i] - fact * temp
            if i < n - 2:
                dl[i] = du[i + 1]
                du[i + 1] = -fact * dl[i]
            du[i] = temp
            temp = x[i]
            x[i] = x[i + 1]
            x[i + 1] = temp - fact * x[i + 1]
    x[n - 1] = x[n - 1] / dd[n - 1]
    x[n - 2] = (x[n - 2] - du[n - 2] * x[n - 1]) / dd[n - 2]
    for i in range(n - 3, -1, -1):
        x[i] = (x[i] - du[i] * x[i + 1] - dl[i] * x[i + 2]) / dd[i]


@_JIT
def _coeff(
    abc, wcfac, gam_fro_fac, efficacy, coeffs
):  # pylint: disable=too-many-locals
    """Fill ``abc`` with one hemisphere's coefficients (``coeff``, efficacy)."""
    cc1, ca0, idz0, idz1, idzl, ca_mid, cc2, idzm = coeffs
    lm = len(ca_mid) + 2
    a, b, c = abc[0], abc[1], abc[2]
    c[0] = 0.0
    c[1] = cc1
    for i in range(2, lm):
        c[i] = cc2[i - 2] - wcfac * idzm[i - 2]
    a[0] = ca0 + wcfac * idz1
    for i in range(1, lm - 1):
        a[i] = ca_mid[i - 1]
    a[lm - 1] = 0.0
    b[0] = 1.0 - c[1] + gam_fro_fac - wcfac * idz0
    for i in range(1, lm - 1):
        b[i] = 1.0 - a[i - 1] - c[i + 1]
    b[lm - 1] = 1.0 - a[lm - 2] + wcfac * idzl
    if efficacy != 1.0:
        c_eff = c[1] * efficacy
        b[0] = b[0] + c[1] - c_eff
        c[1] = c_eff


@_JIT
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def run_years(  # noqa: PLR0913, PLR0915
    state, forc, volc_n, volc_s, rlamda, efficacy, pam, geom, out
):  # pylint: disable=too-many-locals,too-many-statements,invalid-name
    """
    Step the energy budget through ``len(forc)`` years

    Parameters
    ----------
    state : tuple of np.ndarray
        ``(tn, ts, prev)``: ocean temperature profiles and the previous
        year's hemispheric forcing, updated in place
    forc : np.ndarray
        Hemispheric forcing, shape (years, 2)
    volc_n, volc_s : np.ndarray
        Volcanic forcing, shape (years, steps per year or fewer)
    rlamda, efficacy : np.ndarray
        Feedback parameter and ocean heat-uptake efficacy for each
        sub-annual step, shape (years, ldtime)
    pam : tuple of float
        Run-constant parameters (see ``UpwellingDiffusionModel._kernel_args``)
    geom : tuple
        ``(dz, deep_n, deep_s, coeffs)``: layer thicknesses, deep-ocean
        exchange factors and the constant parts of the coefficients
    out : np.ndarray
        Filled with the results named in ``OUT_KEYS``, shape (years, 11)
    """
    tn, ts, prev = state
    W, cpi, beto, threstemp, rlamdo, foan, foas, ebbeta = pam[0:8]
    dt, c1, fnso, year_s, eff0, hcap, ohc_fac, frac700, max700 = pam[8:17]
    dz, deep_n, deep_s, coeffs = geom
    max700 = int(max700)
    lm, ldtime = len(tn), rlamda.shape[1]
    dtyear = 1.0 / ldtime
    fnsa = 1.0
    blm = ebbeta / rlamdo
    c1fac = dt / (c1 * dz[0])
    dtmnl3 = dt * beto / (c1 * dz[lm - 1])
    dtmnl1 = 1.0 - dtmnl3
    dtmsl3 = fnso * dtmnl3
    dtmsl1 = 1.0 - dtmsl3
    dn, ds = np.empty(lm), np.empty(lm)
    abc_n, abc_s = np.empty((3, lm)), np.empty((3, lm))
    for y in range(len(forc)):
        forc_nh, forc_sh = forc[y, 0], forc[y, 1]
        tn_start, ts_start = tn.copy(), ts.copy()
        temp1n = temp1s = 0.0
        tempn = temps = tempn_air = temps_air = tempn_sea = temps_sea = 0.0
        for im in range(ldtime):
            rl = rlamda[y, im]
            fnx = rl + foan * rlamdo + ebbeta
            fsx = rl + foas * rlamdo + ebbeta
            gamn = foan + rl / rlamdo + blm
            gams = foas + rl / rlamdo + fnsa * blm
            den = gams * gamn - fnsa * blm**2
            dtrm2n = (beto + foas * ebbeta / den) * c1fac
            dtrm3n = gams / den * c1fac
            dtrm4n = blm / den * c1fac
            dtrm2s = (fnso * beto + foan * fnsa * ebbeta / den) * c1fac
            dtrm3s = gamn / den * c1fac
            dtrm4s = fnsa * blm / den * c1fac
            gffn = rlamdo * (1.0 - foan * gams / (gams * gamn - blm * blm)) * dt
            gffn = gffn / (c1 * dz[0])
            gffs = rlamdo * (1.0 - foas * gamn / (gamn * gams - blm * blm)) * dt
            gffs = gffs / (c1 * dz[0])
            if threstemp == 0:
                wcn = wcs = W / year_s * dt
            else:
                wcn = W / year_s * (1 - 0.3 * temp1n / threstemp) * dt
                wcs = W / year_s * (1 - 0.3 * temp1s / threstemp) * dt
            _coeff(abc_n, wcn, gffn, efficacy[y, im], coeffs)
            _coeff(abc_s, wcs, gffs, efficacy[y, im], coeffs)
            dtrm1n = 1.0 - cpi * wcn / dz[0] - beto * dt / (c1 * dz[0])
            dtrm1s = 1.0 - cpi * wcs / dz[0] - fnso * beto * dt / (c1 * dz[0])
            dtmnl2 = wcn * cpi / dz[lm - 1]
            dtmsl2 = wcs * cpi / dz[lm - 1]

            vn = volc_n[y, im % volc_n.shape[1]]
            vs = volc_s[y, im % volc_s.shape[1]]
            dqn = (im + 1) * forc_nh * dtyear + (1 - (im + 1) * dtyear) * prev[0] + vn
            dqs = (im + 1) * forc_sh * dtyear + (1 - (im + 1) * dtyear) * prev[1] + vs
            dn[0] = dtrm1n * tn[0] + dtrm2n * ts[0] + dtrm3n * dqn + dtrm4n * dqs
            ds[0] = dtrm1s * ts[0] + dtrm2s * tn[0] + dtrm3s * dqs + dtrm4s * dqn
            for i in range(1, lm - 1):
                dn[i] = tn[i] + deep_n[i - 1] * (ts[i] - tn[i])
                ds[i] = ts[i] + deep_s[i - 1] * (tn[i] - ts[i])
            dn[lm - 1] = dtmnl1 * tn[lm - 1] + dtmnl2 * tn[0] + dtmnl3 * ts[lm - 1]
            ds[lm - 1] = dtmsl1 * ts[lm - 1] + dtmsl2 * ts[0] + dtmsl3 * tn[lm - 1]
            _gtsv(abc_n[0], abc_n[1], abc_n[2], dn)
            _gtsv(abc_s[0], abc_s[1], abc_s[2], ds)
            tn[:], ts[:] = dn, ds
            temp1n, temp1s = tn[0], ts[0]

            tempan = (
                dqn
                + foan * rlamdo * temp1n
                + ebbeta * (dqs + foas * rlamdo * temp1s) / fsx
            )
            tempan = tempan / (fnx - ebbeta**2 / fsx)
            tempas = (
                dqs
                + foas * rlamdo * temp1s
                + ebbeta * (dqn + foan * rlamdo * temp1n) / fnx
            )
            tempas = tempas / (fsx - ebbeta**2 / fnx)
            tempn = tempn + (foan * temp1n + (1.0 - foan) * tempan) / ldtime
            temps = temps + (foas * temp1s + (1.0 - foas) * tempas) / ldtime
            tempn_air = tempn_air + tempan / ldtime
            temps_air = temps_air + tempas / ldtime
            tempn_sea = tempn_sea + temp1n / ldtime
            temps_sea = temps_sea + temp1s / ldtime
        prev[0], prev[1] = forc_nh, forc_sh

        heat = hcap * dz * ((tn - tn_start) * foan + (ts - ts_start) * foas)
        anom = (eff0 - 1.0) * (_sum(heat[1:]) / year_s)
        havtemp = ohc_fac * dz * (tn * foan + ts * foas)
        ohc700 = _sum(havtemp[:max700]) + frac700 * havtemp[max700]
        ribn = forc_nh + _sum(volc_n[y]) / volc_n.shape[1] - rl * tempn - anom
        ribs = forc_sh + _sum(volc_s[y]) / volc_s.shape[1] - rl * temps - anom
        out[y, 0], out[y, 1] = tempn, temps
        out[y, 2], out[y, 3] = tempn_air, temps_air
        out[y, 4], out[y, 5] = tempn_sea, temps_sea
        out[y, 6], out[y, 7], out[y, 8] = ribn, ribs, anom
        out[y, 9], out[y, 10] = ohc700, _sum(havtemp)
