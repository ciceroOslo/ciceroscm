"""5x5 sweep over delta_lambda_pdo x pdo_efficacy_scale on an SSP2-4.5 historical run.

Runs CICEROSCM from 1750 to 2025 with the observed PDO index loaded from
``tests/test-data/pdo_ts_noaa.dat`` (1870-2025, padded with zeros
before 1870). All other UDM parameters held at the integration-test defaults.
Output: a 5x5 grid of dT_glob time series, one subplot per parameter combo.
"""

import os

import matplotlib.pyplot as plt
import numpy as np

from ciceroscm import CICEROSCM

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DATA = os.path.join(REPO_ROOT, "tests", "test-data")

NYSTART = 1750
NYEND = 2025

DELTA_LAMBDA_PDO_VALUES = [-0.5, -0.2, 0.0, 0.2, 0.5]
PDO_EFFICACY_SCALE_VALUES = [-0.3, -0.1, 0.0, 0.1, 0.3]

PAMSET_UDM_BASE = {
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

PAMSET_EMICONC = {
    "qbmb": 0.0,
    "qo3": 0.5,
    "qdirso2": -0.00308,
    "qindso2": -0.97 / 57.052577209999995,
    "qbc": 0.0279,
    "qoc": -0.00433,
    "qh2o_ch4": 0.091915,
    "ref_yr": 2010,
}


def load_pdo_monthly(nystart, nyend):
    """Read pdo_ts_noaa.dat and return a (nyears, 12) array."""
    path = os.path.join(TEST_DATA, "pdo_ts_noaa.dat")
    dates, values = [], []
    with open(path, "r", encoding="utf-8") as fh:
        raw = np.loadtxt(fh, dtype=float, skiprows=1)
        print(raw)
        years = raw[:, 0].astype(int)
        values = raw[:, 1:]

    padded = np.zeros((nyend - nystart + 1, 12), dtype=float)
    mask = (years >= nystart) & (years <= nyend)
    padded[years[mask] - nystart, :] = values[mask]
    padded = np.where(padded == -99.99, 0.0, padded)  # Replace missing value code with 0.0
    return padded


def build_cscm():
    return CICEROSCM(
        {
            "gaspam_file": os.path.join(TEST_DATA, "gases_vupdate_2022_AR6.txt"),
            "nystart": NYSTART,
            "nyend": NYEND,
            "emstart": 1850,
            "concentrations_file": os.path.join(TEST_DATA, "ssp245_conc_RCMIP.txt"),
            "emissions_file": os.path.join(TEST_DATA, "ssp245_em_RCMIP.txt"),
            "nat_ch4_file": os.path.join(TEST_DATA, "natemis_ch4.txt"),
            "nat_n2o_file": os.path.join(TEST_DATA, "natemis_n2o.txt"),
        }
    )


def run_one(pdo_data, delta_lambda_pdo, pdo_efficacy_scale):
    cscm = build_cscm()
    pamset = {
        **PAMSET_UDM_BASE,
        "pdo_index_data": pdo_data,
        "delta_lambda_pdo": delta_lambda_pdo,
        "pdo_efficacy_scale": pdo_efficacy_scale,
    }
    cscm._run(
        {"results_as_dict": True},
        pamset_udm=pamset,
        pamset_emiconc=PAMSET_EMICONC,
    )
    return (
        np.asarray(cscm.results["dT_glob"]),
        np.asarray(cscm.results["RIB_glob"]),
    )


def local_lambda(temp, rib, window=30):
    """OLS slope of RIB on dT in a centered window, following Andrews et al. (2022).

    Returns an array of length ``len(temp)`` with NaN where the centered
    window falls outside the series. Slope is dN/dT (W m^-2 K^-1); the
    Gregory feedback is the negative of this slope.
    """
    half = window // 2
    out = np.full_like(temp, np.nan, dtype=float)
    for t in range(half, len(temp) - (window - half)):
        x = temp[t - half : t - half + window]
        y = rib[t - half : t - half + window]
        xm = x - x.mean()
        denom = (xm * xm).sum()
        if denom > 0:
            out[t] = (xm * (y - y.mean())).sum() / denom
    return out


def main():
    pdo_data = load_pdo_monthly(NYSTART, NYEND)
    years = np.arange(NYSTART, NYEND + 1)

    nrows = len(PDO_EFFICACY_SCALE_VALUES)
    ncols = len(DELTA_LAMBDA_PDO_VALUES)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(3.0 * ncols, 2.2 * nrows), sharex=True, sharey=True
    )

    all_series = np.empty((nrows, ncols, years.size))
    all_rib = np.empty((nrows, ncols, years.size))
    for i, eff in enumerate(PDO_EFFICACY_SCALE_VALUES):
        for j, dlam in enumerate(DELTA_LAMBDA_PDO_VALUES):
            dt_glob, rib = run_one(pdo_data, dlam, eff)
            all_series[i, j] = dt_glob
            all_rib[i, j] = rib

    i0 = PDO_EFFICACY_SCALE_VALUES.index(0.0)
    j0 = DELTA_LAMBDA_PDO_VALUES.index(0.0)
    baseline = all_series[i0, j0]
    anomalies = all_series - baseline
    vmax = np.max(np.abs(anomalies))

    for i, eff in enumerate(PDO_EFFICACY_SCALE_VALUES):
        for j, dlam in enumerate(DELTA_LAMBDA_PDO_VALUES):
            ax = axes[i, j]
            ax.plot(years, anomalies[i, j], lw=0.9, color="C0")
            ax.set_title(
                f"$\\Delta\\lambda_{{PDO}}$={dlam:+.2f}, "
                f"$\\epsilon_{{PDO}}$={eff:+.2f}",
                fontsize=8,
            )
            ax.axhline(0, color="0.7", lw=0.4)
            ax.set_ylim(-vmax * 1.05, vmax * 1.05)
            ax.set_xlim(1870, 2025)
            ax.grid(True, alpha=0.3)

    for ax in axes[-1, :]:
        ax.set_xlabel("Year")
    for ax in axes[:, 0]:
        ax.set_ylabel("dT_glob anomaly (K)")

    fig.suptitle(
        "CICEROSCM SSP2-4.5: PDO-induced dT_glob anomaly vs (0,0) baseline "
        "($\\Delta\\lambda_{PDO}$ cols, $\\epsilon_{PDO}$ rows)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    out_path = os.path.join(REPO_ROOT, "scripts", "pdo_sweep_dT_glob_anomaly.png")
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")

    # ----- second figure: 30-year sliding-window OLS feedback (Andrews et al. 2022) -----
    window = 30
    lam_all = np.empty_like(all_series)
    for i in range(nrows):
        for j in range(ncols):
            lam_all[i, j] = local_lambda(all_series[i, j], all_rib[i, j], window=window)

    lam_anom = lam_all - lam_all[i0, j0]
    finite = np.isfinite(lam_anom)
    vmax_lam = np.nanmax(np.abs(lam_anom[finite]))

    fig2, axes2 = plt.subplots(
        nrows, ncols, figsize=(3.0 * ncols, 2.2 * nrows), sharex=True, sharey=True
    )

    for i, eff in enumerate(PDO_EFFICACY_SCALE_VALUES):
        for j, dlam in enumerate(DELTA_LAMBDA_PDO_VALUES):
            ax = axes2[i, j]
            ax.plot(years, lam_anom[i, j], lw=0.9, color="C3")
            ax.set_title(
                f"$\\Delta\\lambda_{{PDO}}$={dlam:+.2f}, "
                f"$\\epsilon_{{PDO}}$={eff:+.2f}",
                fontsize=8,
            )
            ax.axhline(0, color="0.7", lw=0.4)
            ax.set_xlim(1870, 2025)
            ax.set_ylim(-vmax_lam * 1.05, vmax_lam * 1.05)
            ax.grid(True, alpha=0.3)

    for ax in axes2[-1, :]:
        ax.set_xlabel("Year (centre of window)")
    for ax in axes2[:, 0]:
        ax.set_ylabel(r"$\Delta$(dN/dT) (W m$^{-2}$ K$^{-1}$)")

    fig2.suptitle(
        f"Local feedback anomaly vs (0,0) baseline: {window}-yr centred OLS "
        "slope of RIB_glob vs dT_glob (Andrews 2022 method)",
        fontsize=11,
    )
    fig2.tight_layout(rect=(0, 0, 1, 0.97))

    out_path2 = os.path.join(REPO_ROOT, "scripts", "pdo_sweep_local_lambda_anomaly.png")
    fig2.savefig(out_path2, dpi=150)
    print(f"wrote {out_path2}")


if __name__ == "__main__":
    main()
