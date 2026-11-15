"""Daily NEE 1997-2014: the tower's data as dots, the fitted model as a line with two bands.

Nothing is sampled. Posterior draws are run forward from 1 January 1997 to 31 December 2014,
exactly as in script 06, so the test years 2011-2014 are a true forecast from the fitted years.

    dots         observed daily NEE on the days used (QC at least 0.75 and RANDUNC present, D12)
    line         model median, every day
    dark band    90% band from the parameters only (the spread of the posterior forests)
    light band   90% predictive band: parameters + measurement error (RANDUNC) + model error (sigma_m),
                 the band the test coverage (D61) is computed on. On days without RANDUNC
                 only sigma_m is added.
Top panel: the whole record, fitted years and test years marked. Bottom panel: 2010-2013 enlarged.

Run from the repository root:
    python scripts/17_nee_timeseries.py fit_model_error_MC_cold_dense_ta995_d2000
    python scripts/17_nee_timeseries.py --data-only      (the tower's data alone, no model, no fit needed;
                                                          saves fig_Z_observed_nee_timeseries)

Saves in figures/:  fig_Z_<fit>_nee_timeseries (pdf and png)
Prints the share of used days inside each band, for the fitted and the test years.
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, year_block
from dalec.data_io import MASK_USE, load_record_with_drivers, select_years
from dalec.plotting import BLUE, GREY_LIGHT, INK, INK_SOFT, UNIT_FLUX, apply_style, save_figure

_spec = importlib.util.spec_from_file_location("evaluate", Path(__file__).resolve().parent / "06_evaluate.py")
evaluate = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(evaluate)

ZOOM = ("2010-01-01", "2013-12-31")


def bands(model_nee: np.ndarray, sigma_m: np.ndarray, sigma: np.ndarray, interval: float, seed: int):
    """Median, parameter-only band and full predictive band per day (draws x days in, days out)."""
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    median = np.median(model_nee, axis=0)
    p_lo, p_hi = np.percentile(model_nee, [lo, hi], axis=0)
    obs_sd = np.where(np.isfinite(sigma), sigma, 0.0)
    sd = np.sqrt(obs_sd[None, :] ** 2 + sigma_m[:, None] ** 2)
    replicated = model_nee + sd * np.random.default_rng(seed).standard_normal(model_nee.shape)
    f_lo, f_hi = np.percentile(replicated, [lo, hi], axis=0)
    return median, (p_lo, p_hi), (f_lo, f_hi)


def share_inside(observed: np.ndarray, band: tuple[np.ndarray, np.ndarray], mask: np.ndarray) -> float:
    lo, hi = band
    inside = (observed >= lo) & (observed <= hi)
    return float(inside[mask].mean()) if mask.any() else float("nan")


def draw(ax, dates, observed, use, median, param_band, full_band, test_start, show_labels):
    ax.fill_between(dates, *full_band, color=BLUE, alpha=0.18, lw=0, label="90% band: all uncertainty")
    ax.fill_between(dates, *param_band, color=BLUE, alpha=0.55, lw=0, label="90% band: parameters only")
    ax.plot(dates, median, color=BLUE, lw=0.8, label="model (median)")
    ax.plot(dates[use], observed[use], "o", color=INK, ms=1.4, alpha=0.6, label="tower data (days used)")
    ax.axvspan(pd.Timestamp(f"{test_start}-01-01"), dates[-1], color=GREY_LIGHT, alpha=0.35, lw=0, zorder=0)
    ax.axhline(0, color=INK_SOFT, lw=0.6)
    ax.set_ylabel(f"NEE\n({UNIT_FLUX})")
    ax.grid(False)
    if show_labels:
        top = ax.get_ylim()[1]
        ax.text(dates[0], top, "  fitted years", va="top", color=INK_SOFT, fontsize=10)
        ax.text(pd.Timestamp(f"{test_start}-01-01"), top, "  test years (not fitted)", va="top",
                color=INK_SOFT, fontsize=10)


def figure(dates, observed, use, median, param_band, full_band, test_start):
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(12, 7.2), gridspec_kw={"height_ratios": [1, 1.15]})
    draw(top, dates, observed, use, median, param_band, full_band, test_start, show_labels=True)
    top.set_xlim(dates[0], dates[-1])
    top.xaxis.set_major_locator(mdates.YearLocator(2))
    top.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    zoom = (dates >= pd.Timestamp(ZOOM[0])) & (dates <= pd.Timestamp(ZOOM[1]))
    draw(bottom, dates[zoom], observed[zoom], use[zoom], median[zoom],
         tuple(b[zoom] for b in param_band), tuple(b[zoom] for b in full_band), test_start, show_labels=False)
    bottom.set_xlim(dates[zoom][0], dates[zoom][-1])
    bottom.xaxis.set_major_locator(mdates.YearLocator())
    bottom.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    bottom.set_title(f"{ZOOM[0][:4]} to {ZOOM[1][:4]}, enlarged", loc="left", fontsize=11, color=INK_SOFT)
    handles, labels = top.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncols=4, frameon=False, markerscale=4)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return fig


def data_figure(dates, observed, use, test_start):
    """The tower's daily NEE alone: days used in black, other measured days in light grey."""
    fig, ax = plt.subplots(figsize=(12, 4.2))
    other = np.isfinite(observed) & ~use
    ax.plot(dates[other], observed[other], "o", color=GREY_LIGHT, ms=1.6, label="other measured days")
    ax.plot(dates[use], observed[use], "o", color=INK, ms=1.6, alpha=0.7, label="days used (QC and RANDUNC)")
    ax.axvspan(pd.Timestamp(f"{test_start}-01-01"), dates[-1], color=GREY_LIGHT, alpha=0.35, lw=0, zorder=0)
    ax.axhline(0, color=INK_SOFT, lw=0.6)
    ax.set_xlim(dates[0], dates[-1])
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_ylabel(f"NEE\n({UNIT_FLUX})")
    ax.grid(False)
    top = ax.get_ylim()[1]
    ax.text(dates[0], top, "  fitted years", va="top", color=INK_SOFT, fontsize=10)
    ax.text(pd.Timestamp(f"{test_start}-01-01"), top, "  test years (not fitted)", va="top",
            color=INK_SOFT, fontsize=10)
    fig.legend(loc="lower center", ncols=2, frameon=False, markerscale=4)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fit", nargs="?", help="real-data fit in results/, without .nc")
    parser.add_argument("--data-only", action="store_true", help="draw the tower's data alone, no model")
    args = parser.parse_args()
    config = load_config()
    apply_style()
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 12, "xtick.labelsize": 10.5,
                         "ytick.labelsize": 10.5, "legend.fontsize": 11})
    ev = config["evaluation"]
    n, seed, interval = int(ev["n_draws"]), int(ev["seed"]), float(ev["interval"])
    cal_start, cal_end = year_block(config, "calibration")
    test_start, test_end = year_block(config, "evaluation")
    data = select_years(load_record_with_drivers(config), cal_start, test_end)
    cols = config["columns"]
    observed = data[cols["nee"]].to_numpy(dtype=float)
    sigma = data[cols["sigma"]].to_numpy(dtype=float)
    use = data[MASK_USE].to_numpy(dtype=bool)
    years = data.index.year.to_numpy()

    if args.data_only:
        for name, mask in (("fitted years", use & (years <= cal_end)), ("test years", use & (years >= test_start))):
            print(f"  {name}: {mask.sum()} days used")
        save_figure(data_figure(data.index, observed, use, test_start), "fig_Z_observed_nee_timeseries", config)
        return
    if not args.fit:
        parser.error("give a fit name, or --data-only")

    idata = evaluate.load_fit(args.fit, config)
    passed, numbers = evaluate.computation_checks(idata, config)
    print(f"  D49: divergences {numbers['divergences']}, largest R-hat {numbers['max_rhat']:.4f}, "
          f"smallest ESS {numbers['min_ess']:.0f}  ->  {'PASS' if passed else 'FAIL (figure NOT FOR REPORTING)'}")
    model_nee, sigma_m = evaluate.forward(idata, data, config, args.fit, n, seed)
    median, param_band, full_band = bands(model_nee, sigma_m, sigma, interval, seed + 1)

    for name, mask in (("fitted years", use & (years <= cal_end)), ("test years", use & (years >= test_start))):
        print(f"  {name}: {mask.sum()} days used; inside the parameter band "
              f"{share_inside(observed, param_band, mask):.1%}, inside the full band "
              f"{share_inside(observed, full_band, mask):.1%}")
    fig = figure(data.index, observed, use, median, param_band, full_band, test_start)
    save_figure(fig, f"fig_Z_{args.fit}_nee_timeseries", config)


if __name__ == "__main__":
    main()
