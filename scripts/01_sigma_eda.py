"""Exploratory figures for the observation uncertainty (sigma = RANDUNC).

Makes five figures for the calibration block (1997-2010):

    A  fig_A_sigma_over_time     RANDUNC for every day; where the gaps are
    B  fig_B_sigma_vs_nee        RANDUNC against |NEE|, with a straight trend line
    C  fig_C_sigma_by_month      box plots of RANDUNC for each calendar month
    D  fig_D_sigma_distribution  histogram and smoothed density of RANDUNC
    E  fig_E_sigma_annual_cycle  typical RANDUNC for each day of the year

Which days each figure uses (DECISIONS.md, D15):
    A uses every day of the block, because it shows what is excluded and missing.
    B to E use assimilable days only, because those are what the likelihood sees.

Run from the repository root:
    python scripts/01_sigma_eda.py
"""

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from dalec.config import load_config, year_block
from dalec.data_io import MASK_USE, load_screened_block, screening_report
from dalec.plotting import (
    BLUE, GREY, GREY_LIGHT, INK, MONTHS, UNIT_FLUX,
    apply_style, climatology, heading, save_figure, set_year_axis,
)


def missing_runs(is_missing: pd.Series) -> pd.DataFrame:
    """Find stretches of consecutive True days. Returns start, end and length."""
    flags = is_missing.to_numpy(dtype=bool)
    edges = np.diff(np.concatenate([[0], flags.astype(int), [0]]))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1
    dates = is_missing.index
    return pd.DataFrame({
        "start": dates[starts].date,
        "end": dates[ends].date,
        "days": ends - starts + 1,
    })


# ---------------------------------------------------------------------------
# Figure A: sigma over time (every day)
# ---------------------------------------------------------------------------

def figure_a(data: pd.DataFrame, config: dict[str, Any]) -> None:
    """RANDUNC for every day, split into assimilable, excluded, and missing."""
    sigma = data[config["columns"]["sigma"]]
    use = data[MASK_USE]
    has_sigma = sigma.notna() & (sigma > 0)
    excluded = has_sigma & ~use

    heading("Figure A: RANDUNC over time")
    print(f"  assimilable days               {int(use.sum()):>5}")
    print(f"  RANDUNC present but excluded   {int(excluded.sum()):>5}")
    print(f"  RANDUNC missing or zero        {int((~has_sigma).sum()):>5}")
    runs = missing_runs(~has_sigma)
    long_runs = runs[runs["days"] >= 7].sort_values("days", ascending=False)
    print(f"\n  Gaps of 7 days or more ({len(long_runs)} gaps, "
          f"{int(long_runs['days'].sum())} of {int(runs['days'].sum())} missing days):")
    print(long_runs.to_string(index=False) if len(long_runs) else "  none")

    fig, ax = plt.subplots(figsize=(10, 3.4))
    for _, run in runs.iterrows():
        ax.axvspan(pd.Timestamp(run["start"]), pd.Timestamp(run["end"]) + pd.Timedelta(days=1),
                   color=GREY_LIGHT, lw=0, zorder=0)
    ax.scatter(sigma.index[excluded], sigma[excluded], s=4, color=GREY,
               lw=0, label="RANDUNC present, excluded (QC below threshold)", zorder=2)
    ax.scatter(sigma.index[use], sigma[use], s=4, color=BLUE,
               lw=0, label="assimilable", zorder=3)
    ax.fill_between([], [], color=GREY_LIGHT, label="RANDUNC missing")  # legend entry only

    start, end = year_block(config, "calibration")
    ax.set_title(f"Daily random uncertainty of NEE (RANDUNC), FI-Hyy {start}-{end}")
    ax.set_ylabel(f"RANDUNC ({UNIT_FLUX})")
    ax.set_xlim(pd.Timestamp(start, 1, 1), pd.Timestamp(end, 12, 31))
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncols=3, markerscale=3)
    save_figure(fig, "fig_A_sigma_over_time", config)


# ---------------------------------------------------------------------------
# Figure B: sigma against |NEE| (assimilable days)
# ---------------------------------------------------------------------------

def figure_b(sigma: pd.Series, abs_nee: pd.Series, config: dict[str, Any]) -> None:
    """RANDUNC against flux size, one colour, with a straight trend line."""
    slope, intercept = np.polyfit(abs_nee, sigma, deg=1)

    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.scatter(abs_nee, sigma, s=6, color=BLUE, alpha=0.35, lw=0, label="assimilable days")
    xs = np.linspace(0, abs_nee.max(), 50)
    ax.plot(xs, intercept + slope * xs, color=INK, lw=1.5, label="straight trend line")

    ax.set_title("RANDUNC against the size of the daily flux")
    ax.set_xlabel(f"|NEE| ({UNIT_FLUX})")
    ax.set_ylabel(f"RANDUNC ({UNIT_FLUX})")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left")
    save_figure(fig, "fig_B_sigma_vs_nee", config)


# ---------------------------------------------------------------------------
# Figure C: sigma by calendar month (assimilable days)
# ---------------------------------------------------------------------------

def figure_c(sigma: pd.Series, config: dict[str, Any]) -> None:
    """Box plots of RANDUNC for each calendar month, all years pooled."""
    month = sigma.index.month
    groups = [sigma[month == m].to_numpy() for m in range(1, 13)]

    summary = pd.DataFrame({
        "days": [len(g) for g in groups],
        "median": [np.median(g) for g in groups],
        "q25": [np.percentile(g, 25) for g in groups],
        "q75": [np.percentile(g, 75) for g in groups],
    }, index=MONTHS)
    heading("Figure C: RANDUNC by calendar month (assimilable days)")
    print(summary.round(3).to_string())

    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.boxplot(
        groups,
        tick_labels=MONTHS,
        widths=0.6,
        patch_artist=True,
        boxprops={"facecolor": BLUE, "edgecolor": BLUE, "alpha": 0.85},
        medianprops={"color": "white", "linewidth": 1.5},
        whiskerprops={"color": BLUE},
        capprops={"color": BLUE},
        flierprops={"marker": "o", "markersize": 2, "markerfacecolor": GREY,
                    "markeredgewidth": 0},
    )
    start, end = year_block(config, "calibration")
    ax.set_title(f"RANDUNC by calendar month, all years {start}-{end}")
    ax.set_ylabel(f"RANDUNC ({UNIT_FLUX})")
    ax.set_ylim(bottom=0)
    ax.grid(axis="x", visible=False)
    save_figure(fig, "fig_C_sigma_by_month", config)


# ---------------------------------------------------------------------------
# Figure D: distribution of sigma (assimilable days)
# ---------------------------------------------------------------------------

def figure_d(sigma: pd.Series, config: dict[str, Any]) -> None:
    """Histogram of RANDUNC with a smoothed density curve on top."""
    values = sigma.to_numpy()
    median = float(np.median(values))

    heading("Figure D: distribution of RANDUNC (assimilable days)")
    print(f"  median {median:.3f}   mean {values.mean():.3f}   "
          f"5th-95th percentile {np.percentile(values, 5):.3f} to {np.percentile(values, 95):.3f}")

    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.hist(values, bins=60, density=True, color=BLUE, alpha=0.45,
            edgecolor="white", linewidth=0.3, label="histogram")
    # Draw the curve only between the smallest and largest observed value:
    # smoothing would otherwise show density below the minimum, where no data exist.
    xs = np.linspace(values.min(), values.max(), 400)
    ax.plot(xs, stats.gaussian_kde(values)(xs), color=INK, lw=1.5, label="smoothed density")
    ax.axvline(median, color=INK, lw=1, ls=":")
    ax.text(median, 0.98, f" median {median:.2f}", transform=ax.get_xaxis_transform(),
            va="top", ha="left", fontsize=8)

    ax.set_title("How often each RANDUNC value occurs")
    ax.set_xlabel(f"RANDUNC ({UNIT_FLUX})")
    ax.set_ylabel("density")
    ax.set_xlim(left=0)
    ax.legend(loc="upper right")
    save_figure(fig, "fig_D_sigma_distribution", config)


# ---------------------------------------------------------------------------
# Figure E: mean annual cycle (assimilable days)
# ---------------------------------------------------------------------------

def figure_e(sigma: pd.Series, config: dict[str, Any]) -> None:
    """Typical RANDUNC for each day of the year, across all calibration years.

    Uses the median (the middle value), because RANDUNC is skewed: a few very
    large values would pull a mean upwards. The band shows the middle half of
    the values (25th to 75th percentile) for each day.
    """
    cycle = climatology(sigma)

    heading("Figure E: annual cycle of RANDUNC (assimilable days)")
    print(f"  years contributing per day: min {int(cycle['n'].min())}, max {int(cycle['n'].max())}")
    low, high = cycle["median"].idxmin(), cycle["median"].idxmax()
    print(f"  lowest median  {cycle['median'].min():.3f} on day {int(low)}")
    print(f"  highest median {cycle['median'].max():.3f} on day {int(high)}")

    fig, ax = plt.subplots(figsize=(7, 3.8))
    ax.fill_between(cycle.index, cycle["q25"], cycle["q75"], color=BLUE, alpha=0.25,
                    lw=0, label="middle half of values (25th to 75th percentile)")
    ax.plot(cycle.index, cycle["median"], color=BLUE, lw=1.5, label="median")

    set_year_axis(ax)
    ax.set_ylim(bottom=0)
    start, end = year_block(config, "calibration")
    ax.set_title(f"Typical RANDUNC through the year, {start}-{end}")
    ax.set_ylabel(f"RANDUNC ({UNIT_FLUX})")
    ax.legend(loc="upper left")
    save_figure(fig, "fig_E_sigma_annual_cycle", config)


# ---------------------------------------------------------------------------

def main() -> None:
    config = load_config()
    apply_style()

    data = load_screened_block(config, block="calibration")
    cols = config["columns"]
    use = data[MASK_USE]
    sigma = data.loc[use, cols["sigma"]]          # assimilable days only
    abs_nee = data.loc[use, cols["nee"]].abs()

    heading("Screening report (calibration block)")
    print(screening_report(data).to_string())

    figure_a(data, config)
    figure_b(sigma, abs_nee, config)
    figure_c(sigma, config)
    figure_d(sigma, config)
    figure_e(sigma, config)


if __name__ == "__main__":
    main()