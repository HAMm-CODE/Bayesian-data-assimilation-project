"""Shared figure style and helpers for all EDA scripts.

Kept in one place so every figure in the thesis looks the same and uses the
same day-of-year and climatology rules.
"""

from typing import Any

import matplotlib

matplotlib.use("Agg")  # draw to files only, never open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import resolve_path

# Units, written for matplotlib labels.
UNIT_FLUX = "g C m$^{-2}$ d$^{-1}$"

# Colours. Categorical ones were checked for colour-blind safety as a set.
BLUE = "#2a78d6"
ORANGE = "#eb6834"
VIOLET = "#4a3aa7"
INK = "#0b0b0b"          # lines drawn on top of data
INK_SOFT = "#52514e"     # secondary lines and notes
GREY = "#a9a8a2"         # neutral: present but not the focus
GREY_LIGHT = "#e6e5e0"   # bands and backgrounds

MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
MONTH_STARTS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]  # on the 365-day axis

# W m-2 (daily mean) to MJ m-2 d-1 (daily total): 86 400 s per day / 1 000 000 J per MJ.
W_M2_TO_MJ_M2_DAY = 0.0864


def apply_style() -> None:
    """Quiet, print-friendly defaults for every figure."""
    plt.rcParams.update({
        "figure.dpi": 100,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,  # grid behind the data, never on top
        "grid.color": GREY_LIGHT,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "legend.fontsize": 8,
    })


def save_figure(fig: plt.Figure, name: str, config: dict[str, Any]) -> None:
    """Save one figure in every format listed in the config."""
    out_dir = resolve_path(config["paths"]["figures"])
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in config["eda"]["figure_formats"]:
        path = out_dir / f"{name}.{fmt}"
        fig.savefig(path, dpi=config["eda"]["figure_dpi"], bbox_inches="tight")
        print(f"  saved {path}")
    plt.close(fig)


def heading(text: str) -> None:
    """Print a visible section title in the console output."""
    print(f"\n{'=' * 70}\n{text}\n{'=' * 70}")


def day_of_year_365(index: pd.DatetimeIndex) -> pd.Series:
    """Day of year on a 365-day calendar: 29 February is dropped (NaN).

    In leap years every day after 28 February moves back by one, so that
    1 March is always day 60 and 31 December is always day 365.
    """
    doy = pd.Series(index.dayofyear, index=index, dtype=float)
    after_feb = index.is_leap_year & (index.month > 2)
    doy[after_feb] -= 1
    doy[(index.month == 2) & (index.day == 29)] = np.nan
    return doy


def climatology(values: pd.Series) -> pd.DataFrame:
    """Median and middle half (25th to 75th percentile) for each day of the year.

    Missing values are ignored. Returns one row per day 1..365.
    """
    doy = day_of_year_365(values.index)
    keep = doy.notna() & values.notna()
    by_day = values[keep].groupby(doy[keep])
    return pd.DataFrame({
        "median": by_day.median(),
        "q25": by_day.quantile(0.25),
        "q75": by_day.quantile(0.75),
        "n": by_day.size(),
    })


def smooth_cyclic(values: pd.Series, window: int) -> pd.Series:
    """Centred moving average over days 1..365 that wraps December into January.

    Wrapping matters at the year boundary: 31 December is averaged with early
    January, not cut off.
    """
    pad = window
    padded = pd.concat([values.iloc[-pad:], values, values.iloc[:pad]])
    smoothed = padded.rolling(window, center=True, min_periods=1).mean()
    return smoothed.iloc[pad:-pad]


def set_year_axis(ax: plt.Axes) -> None:
    """Month names on a 365-day x-axis."""
    ax.set_xticks(MONTH_STARTS, MONTHS)
    ax.set_xlim(1, 365)