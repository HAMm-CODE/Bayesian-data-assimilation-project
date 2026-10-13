"""Exploratory figures for NEE and the weather drivers.

Makes five figures:

    F   fig_F_data_availability    every day of 1996-2014: assimilable or why not
    G   fig_G_nee_seasonal_cycle   NEE through the year, one line per year coloured
                                   by block, with the spring onset marked
    H   fig_H_model_drivers        the drivers DALEC2 reads, through the year
    I   fig_I_nee_vs_drivers       NEE against sunlight, air and soil temperature,
                                   spring days highlighted
    K   fig_K_nee_autocorrelation  how strongly one day's NEE anomaly predicts the next

Which days each figure uses (DECISIONS.md, D16):
    F  every day of the whole record.
    G  QC-passing days, 1997-2014 (sigma is not needed to describe NEE's shape,
       and requiring it would remove 2014 entirely).
    H  every day of the calibration block (drivers are gap-filled).
    I  QC-passing days of the calibration block.
    K  assimilable days of the calibration block (what the likelihood sees).

Run from the repository root:
    python scripts/02_data_eda.py
"""

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from dalec.config import load_config, resolve_path, year_block
from dalec.data_io import (
    MASK_NEE, MASK_QC, MASK_SIGMA, MASK_USE, add_masks, attach_temperature_extremes,
    daily_temperature_extremes, load_fluxnet_dd,
)
from dalec.plotting import (
    BLUE, GREY, GREY_LIGHT, INK, INK_SOFT, ORANGE, UNIT_FLUX, VIOLET,
    W_M2_TO_MJ_M2_DAY, apply_style, climatology, day_of_year_365, heading, save_figure,
    set_year_axis, smooth_cyclic,
)

UNIT_SW = "MJ m$^{-2}$ d$^{-1}$"

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_whole_record(config: dict[str, Any]) -> pd.DataFrame:
    """Every day in the file, with the likelihood masks and daily max/min temperature added.

    Uses the same building blocks as load_screened_block, without selecting
    a year block, because Figure F shows the whole record. The daily highest and
    lowest temperature come from the half-hourly file (DECISIONS.md, D18).
    """
    cols = config["columns"]
    drv = config["drivers"]
    ext = drv["temperature_extremes"]
    missing = float(config["data"]["missing_value"])
    data = load_fluxnet_dd(
        resolve_path(config["paths"]["fluxnet_dd"]),
        timestamp_column=cols["timestamp"],
        required_columns=[cols["nee"], cols["sigma"], cols["qc"], *drv["model"].values(),
                          ext["daily_check_column"]],
        missing_value=missing,
    )
    extremes = daily_temperature_extremes(
        resolve_path(config["paths"]["fluxnet_hh"]),
        timestamp_column=ext["hh_timestamp"],
        temperature_column=ext["hh_column"],
        half_hours_per_day=int(ext["half_hours_per_day"]),
        missing_value=missing,
    )
    data = attach_temperature_extremes(
        data,
        extremes,
        mean_column=ext["daily_check_column"],
        t_max_name=ext["t_max"],
        t_min_name=ext["t_min"],
        t_mid_name=ext["t_mid"],
        half_hours_per_day=int(ext["half_hours_per_day"]),
        mean_tolerance=float(ext["mean_tolerance"]),
    )
    return add_masks(
        data,
        nee_column=cols["nee"],
        sigma_column=cols["sigma"],
        qc_column=cols["qc"],
        qc_threshold=float(config["data"]["qc_threshold"]),
    )


def in_years(data: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Rows from the inclusive year range [start, end]."""
    years = data.index.year
    return data.loc[(years >= start) & (years <= end)]


def soil_layers(data: pd.DataFrame, prefix: str) -> list[str]:
    """Soil columns for one prefix, shallowest first; QC columns left out."""
    found = [c for c in data.columns if c.startswith(prefix) and not c.endswith("_QC")]
    return sorted(found, key=lambda c: int(c.removeprefix(prefix)))


def spring_onset(nee: pd.Series, qc_pass: pd.Series, run_days: int) -> float:
    """Day of year (365-day axis) on which spring uptake starts in one year.

    Onset = first day of `run_days` days in a row with net uptake (NEE < 0)
    on QC-passing days. A day that fails QC breaks the run, so a gap can
    delay the onset but can never create one. Returns NaN if no such run.
    """
    uptake = (nee < 0) & qc_pass
    full_run = uptake.astype(int).rolling(run_days).sum() == run_days
    if not full_run.any():
        return np.nan
    first_end = full_run.idxmax()                      # date the first full run ends
    start = nee.index[nee.index.get_loc(first_end) - (run_days - 1)]
    doy = day_of_year_365(pd.DatetimeIndex([start])).iloc[0]
    return float(doy) if np.isfinite(doy) else 59.0    # 29 Feb counted as 28 Feb


# ---------------------------------------------------------------------------
# Figure F: data availability, whole record
# ---------------------------------------------------------------------------

def figure_f(data: pd.DataFrame, config: dict[str, Any]) -> None:
    """One row per year, one column per day: assimilable, failed QC, or no RANDUNC."""
    category = pd.Series(np.nan, index=data.index)
    category[data[MASK_USE]] = 0                                    # assimilable
    category[data[MASK_NEE] & ~data[MASK_QC]] = 1                   # QC below threshold
    category[data[MASK_NEE] & data[MASK_QC] & ~data[MASK_SIGMA]] = 2  # QC fine, no RANDUNC

    doy = day_of_year_365(data.index)
    years = list(range(data.index.year.min(), data.index.year.max() + 1))
    grid = np.full((len(years), 365), np.nan)
    keep = doy.notna()
    rows = data.index.year[keep] - years[0]
    grid[rows, doy[keep].astype(int) - 1] = category[keep].to_numpy()

    labels = ["assimilable", "QC below threshold", "QC passes, RANDUNC missing"]
    table = pd.DataFrame({lab: (category == i).groupby(data.index.year).sum()
                          for i, lab in enumerate(labels)})
    heading("Figure F: days per year by outcome (whole record)")
    print(table.astype(int).to_string())

    fig, ax = plt.subplots(figsize=(10, 4.6))
    cmap = ListedColormap([BLUE, GREY, ORANGE])
    cmap.set_bad("white")
    ax.imshow(np.ma.masked_invalid(grid), cmap=cmap, vmin=-0.5, vmax=2.5,
              aspect="auto", interpolation="nearest",
              extent=(0.5, 365.5, len(years) - 0.5, -0.5))

    cal_start, cal_end = year_block(config, "calibration")
    ev_start, ev_end = year_block(config, "evaluation")
    blocks = [("calibration", cal_start, cal_end), ("evaluation", ev_start, ev_end)]
    if years[0] < cal_start:
        blocks.insert(0, ("excluded", years[0], cal_start - 1))
    for name, start, end in blocks:
        top, bottom = start - years[0] - 0.5, end - years[0] + 0.5
        ax.axhline(bottom, color=INK, lw=1)
        span = f"{start}" if start == end else f"{start}-{end}"
        ax.text(372, (top + bottom) / 2, f"{name}\n{span}", va="center", ha="left",
                fontsize=8, clip_on=False)

    ax.set_yticks(range(len(years)), years)
    set_year_axis(ax)
    ax.set_xlim(0.5, 365.5)
    ax.grid(False)
    ax.set_title("What happened to each day of the record")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in (BLUE, GREY, ORANGE)]
    ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncols=3)
    save_figure(fig, "fig_F_data_availability", config)


# ---------------------------------------------------------------------------
# Figure G: NEE seasonal cycle, one line per year, coloured by block
# ---------------------------------------------------------------------------

def figure_g(data: pd.DataFrame, config: dict[str, Any]) -> None:
    """Each year's NEE cycle (QC-passing days), the calibration median, and spring onset."""
    nee_col = config["columns"]["nee"]
    window = int(config["eda"]["smoothing_days"])
    run_days = int(config["eda"]["onset_run_days"])
    cal_start, cal_end = year_block(config, "calibration")
    ev_start, ev_end = year_block(config, "evaluation")

    block = in_years(data, cal_start, ev_end)
    qc_pass = block[MASK_NEE] & block[MASK_QC]
    nee_qc = block[nee_col].where(qc_pass)

    onsets = {}
    fig, ax = plt.subplots(figsize=(8, 4.4))
    for year in range(cal_start, ev_end + 1):
        sel = block.index.year == year
        onsets[year] = spring_onset(block.loc[sel, nee_col], qc_pass[sel], run_days)
        doy = day_of_year_365(block.index[sel])
        line = nee_qc[sel].groupby(doy).first().reindex(range(1, 366))
        line = line.rolling(window, center=True, min_periods=window // 2 + 1).mean()
        is_eval = year >= ev_start
        ax.plot(line.index, line, lw=0.8, alpha=0.8 if is_eval else 0.45,
                color=ORANGE if is_eval else BLUE, zorder=3 if is_eval else 2)

    cal_cycle = climatology(in_years(nee_qc.to_frame(), cal_start, cal_end)[nee_col])
    cal_median = smooth_cyclic(cal_cycle["median"].reindex(range(1, 366)), window)
    ax.plot(cal_median.index, cal_median, color=INK, lw=2, zorder=4)

    onset = pd.Series(onsets, name="onset_doy")
    cal_onsets = onset.loc[cal_start:cal_end].dropna()
    if len(cal_onsets):
        ax.axvspan(cal_onsets.min(), cal_onsets.max(), color=GREY_LIGHT, lw=0, zorder=0)
        ax.axvline(cal_onsets.median(), color=INK_SOFT, lw=1, ls="--", zorder=1)
    ax.axhline(0, color=INK_SOFT, lw=0.8)

    heading("Figure G: spring onset per year (first day of "
            f"{run_days} days in a row with NEE < 0, QC-passing days)")
    print(onset.to_frame().assign(block=lambda t: np.where(t.index >= ev_start,
                                                          "evaluation", "calibration"))
          .to_string())
    if len(cal_onsets):
        print(f"\n  calibration years: median day {cal_onsets.median():.0f}, "
              f"range {cal_onsets.min():.0f} to {cal_onsets.max():.0f}")

    legend = [
        plt.Line2D([], [], color=BLUE, lw=1, alpha=0.6, label=f"calibration years ({cal_start}-{cal_end})"),
        plt.Line2D([], [], color=ORANGE, lw=1, label=f"evaluation years ({ev_start}-{ev_end})"),
        plt.Line2D([], [], color=INK, lw=2, label="calibration median"),
        plt.Line2D([], [], color=INK_SOFT, lw=1, ls="--", label="median spring onset"),
        plt.Rectangle((0, 0), 1, 1, color=GREY_LIGHT, label="range of spring onset"),
    ]
    set_year_axis(ax)
    ax.set_ylabel(f"NEE ({UNIT_FLUX})\nnegative = forest takes up carbon")
    ax.set_title(f"Net ecosystem exchange through the year ({window}-day moving average)")
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncols=3)
    save_figure(fig, "fig_G_nee_seasonal_cycle", config)


# ---------------------------------------------------------------------------
# Figure H: model drivers
# ---------------------------------------------------------------------------

def figure_h(cal: pd.DataFrame, config: dict[str, Any]) -> None:
    """Temperature, sunlight and CO2: the weather inputs DALEC2 reads."""
    drv = config["drivers"]["model"]
    ext = config["drivers"]["temperature_extremes"]
    t_max, t_min, t_mid = ext["t_max"], ext["t_min"], ext["t_mid"]
    sw = cal[drv["sw_in"]] * W_M2_TO_MJ_M2_DAY
    sw_pot = cal[config["eda"]["figure_only_columns"]["sw_in_potential"]] * W_M2_TO_MJ_M2_DAY

    heading("Figure H: model drivers (calibration block)")
    for name in [*drv.values(), t_max, t_min, t_mid]:
        print(f"  {name:12s} missing days: {int(cal[name].isna().sum())}")
    full = int(ext["half_hours_per_day"])
    short = cal["t_half_hours_valid"] < full
    print(f"  days with fewer than {full} valid half-hours: {int(short.sum())}")
    print(f"  mean daily range (max - min): {(cal[t_max] - cal[t_min]).mean():.1f} degC")
    diff = cal[t_mid] - cal[ext["daily_check_column"]]
    print(f"  midpoint minus 24-hour mean ({ext['daily_check_column']}): "
          f"average {diff.mean():+.2f} degC, middle 90% {diff.quantile(0.05):+.2f} to "
          f"{diff.quantile(0.95):+.2f} degC")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))

    ax = axes[0]
    for col, colour, label in ((t_max, ORANGE, "daily highest"),
                               (t_mid, BLUE, "midpoint (drives respiration)"),
                               (t_min, VIOLET, "daily lowest")):
        ax.plot(climatology(cal[col])["median"], color=colour, lw=1.5, label=label)
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    set_year_axis(ax)
    ax.set_ylabel("air temperature (°C)")
    ax.set_title("Air temperature (median per day)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncols=1)

    ax = axes[1]
    clim = climatology(sw)
    ax.fill_between(clim.index, clim["q25"], clim["q75"], color=BLUE, alpha=0.25, lw=0,
                    label="middle half of values")
    ax.plot(clim["median"], color=BLUE, lw=1.5, label="median")
    ax.plot(climatology(sw_pot)["median"], color=INK_SOFT, lw=1, ls="--",
            label="with no atmosphere (potential)")
    set_year_axis(ax)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(f"incoming sunlight ({UNIT_SW})")
    ax.set_title("Incoming sunlight")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncols=1)

    ax = axes[2]
    co2 = cal[drv["co2"]]
    ax.plot(co2.index, co2, color=BLUE, lw=0.6)
    ax.set_ylabel("CO$_2$ (ppm)")
    ax.set_title("Atmospheric CO$_2$ over time")

    fig.tight_layout()
    save_figure(fig, "fig_H_model_drivers", config)


# ---------------------------------------------------------------------------
# Figure I: NEE against sunlight, air and soil temperature
# ---------------------------------------------------------------------------

def figure_i(cal: pd.DataFrame, config: dict[str, Any]) -> None:
    """Scatter of NEE against three drivers, spring days highlighted."""
    nee_col = config["columns"]["nee"]
    drv = config["drivers"]["model"]
    t_mid = config["drivers"]["temperature_extremes"]["t_mid"]
    qc_pass = cal[MASK_NEE] & cal[MASK_QC]
    days = cal.loc[qc_pass]
    spring = days.index.month.isin(config["eda"]["spring_months"])

    x_vars = [
        (days[drv["sw_in"]] * W_M2_TO_MJ_M2_DAY, f"incoming sunlight ({UNIT_SW})"),
        (days[t_mid], "air temperature, midpoint of max and min (°C)"),
    ]
    soil = soil_layers(cal, config["eda"]["figure_only_columns"]["soil_temperature_prefix"])
    if soil:
        x_vars.append((days[soil[0]], "soil temperature, shallowest layer (°C)"))

    heading("Figure I: NEE against drivers (QC-passing days, calibration block)")
    print(f"  spring days {int(spring.sum())}, other days {int((~spring).sum())}")

    fig, axes = plt.subplots(1, len(x_vars), figsize=(4 * len(x_vars), 3.8), sharey=True)
    for ax, (x, label) in zip(np.atleast_1d(axes), x_vars):
        ax.scatter(x[~spring], days.loc[~spring, nee_col], s=4, color=GREY, alpha=0.4, lw=0,
                   label="other months")
        ax.scatter(x[spring], days.loc[spring, nee_col], s=6, color=BLUE, alpha=0.6, lw=0,
                   label="spring (Mar-May)")
        ax.axhline(0, color=INK_SOFT, lw=0.8)
        ax.set_xlabel(label)
    first = np.atleast_1d(axes)[0]
    first.set_ylabel(f"NEE ({UNIT_FLUX})\nnegative = forest takes up carbon")
    first.legend(loc="lower left", markerscale=3)
    fig.suptitle("Daily NEE against its drivers, spring highlighted", fontsize=10)
    fig.tight_layout()
    save_figure(fig, "fig_I_nee_vs_drivers", config)


# ---------------------------------------------------------------------------
# Figure K: day-to-day dependence of NEE anomalies
# ---------------------------------------------------------------------------

def figure_k(cal: pd.DataFrame, config: dict[str, Any]) -> None:
    """Autocorrelation of NEE anomalies on assimilable days.

    Anomaly = a day's NEE minus the typical NEE for that day of the year
    (calibration median, smoothed). Removing the seasonal cycle first stops
    'summer follows summer' from looking like day-to-day dependence.
    """
    nee_col = config["columns"]["nee"]
    window = int(config["eda"]["smoothing_days"])
    max_lag = int(config["eda"]["max_lag_days"])

    nee = cal[nee_col].where(cal[MASK_USE])
    typical = smooth_cyclic(climatology(nee)["median"].reindex(range(1, 366)), window)
    doy = day_of_year_365(cal.index)
    anomaly = nee - doy.map(typical)

    lags = np.arange(1, max_lag + 1)
    acf = np.array([anomaly.autocorr(lag=int(k)) for k in lags])
    pairs = int((anomaly.notna() & anomaly.shift(1).notna()).sum())
    band = 1.96 / np.sqrt(pairs)

    heading("Figure K: autocorrelation of NEE anomalies (assimilable days)")
    for k in (1, 2, 7, 14, max_lag):
        if k <= max_lag:
            print(f"  gap of {k:>2} days: {acf[k - 1]:+.2f}")
    print(f"  range expected if days were independent: +/-{band:.2f}")

    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.axhspan(-band, band, color=GREY_LIGHT, lw=0, label="range expected if days were independent")
    ax.vlines(lags, 0, acf, color=BLUE, lw=2)
    ax.plot(lags, acf, "o", color=BLUE, ms=4)
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    ax.set_xlim(0.5, max_lag + 0.5)
    ax.set_xlabel("gap between the two days (days)")
    ax.set_ylabel("correlation of anomalies")
    ax.set_title("Does a day's unusual NEE carry over to the following days?")
    ax.legend(loc="upper right")
    save_figure(fig, "fig_K_nee_autocorrelation", config)


# ---------------------------------------------------------------------------

def main() -> None:
    config = load_config()
    apply_style()

    record = load_whole_record(config)
    cal = in_years(record, *year_block(config, "calibration"))

    figure_f(record, config)
    figure_g(record, config)
    figure_h(cal, config)
    figure_i(cal, config)
    figure_k(cal, config)


if __name__ == "__main__":
    main()