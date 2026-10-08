"""Load and screen the FLUXNET2015 FULLSET daily (DD) file.

Rules followed throughout:
1. Missing values are made explicit: FLUXNET's -9999 becomes NaN first.
2. No day is ever dropped; screening adds True/False columns instead.
3. Every filter is counted, so thesis numbers trace back to code.

Units: the daily file reports NEE and RANDUNC in g C m-2 d-1 already.
"""

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path, year_block

# Names of the True/False columns added by add_masks().
MASK_NEE = "nee_valid"      # NEE is present
MASK_QC = "qc_pass"         # QC fraction >= threshold
MASK_SIGMA = "sigma_valid"  # RANDUNC is present and > 0
MASK_USE = "assimilable"    # all three: the day enters the likelihood


def load_fluxnet_dd(
    path: str | Path,
    *,
    timestamp_column: str,
    required_columns: list[str],
    missing_value: float = -9999.0,
) -> pd.DataFrame:
    """Read the daily FLUXNET csv into a date-indexed table with NaN for missing."""
    csv_path = Path(path)
    if not csv_path.is_file():
        raise FileNotFoundError(
            f"FLUXNET file not found: {csv_path}\n"
            "Expected the FULLSET daily (DD) csv in data/raw/."
        )

    frame = pd.read_csv(csv_path)

    absent = [c for c in [timestamp_column, *required_columns] if c not in frame.columns]
    if absent:
        raise KeyError(
            f"{csv_path.name} is missing columns {absent}. "
            "Is this the FULLSET (not SUBSET) daily file?"
        )

    non_numeric = [c for c in frame.columns if not pd.api.types.is_numeric_dtype(frame[c])]
    if non_numeric:
        raise TypeError(f"expected only numbers, but these columns hold text: {non_numeric}")

    dates = pd.to_datetime(
        frame[timestamp_column].astype("int64").astype(str), format="%Y%m%d"
    )

    data = frame.drop(columns=[timestamp_column])
    is_missing = np.isclose(data.to_numpy(dtype=float), missing_value, rtol=0.0, atol=1e-6)
    data = data.mask(is_missing)

    data.index = pd.DatetimeIndex(dates, name="date")
    data = data.sort_index()
    _check_daily_contiguous(data.index)
    return data


def _check_daily_contiguous(index: pd.DatetimeIndex) -> None:
    """Stop if any date is duplicated or missing from the calendar.

    The forward model steps one day per row, so a missing or repeated date
    would shift every later day.
    """
    if index.has_duplicates:
        dupes = index[index.duplicated()].unique()[:5]
        shown = ", ".join(str(d) for d in dupes.date)
        raise ValueError(f"the same date appears twice, e.g. {shown}")

    expected = pd.date_range(index.min(), index.max(), freq="D")
    if len(expected) != len(index):
        absent = expected.difference(index)[:5]
        shown = ", ".join(str(d) for d in absent.date)
        raise ValueError(f"dates missing from the calendar, e.g. {shown}")


def select_years(data: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Keep the inclusive year range [start, end], which must be fully covered."""
    years = data.index.year
    selected = data.loc[(years >= start) & (years <= end)]
    if selected.empty:
        raise ValueError(f"no days between {start} and {end}")

    first, last = selected.index.min(), selected.index.max()
    if first != pd.Timestamp(start, 1, 1) or last != pd.Timestamp(end, 12, 31):
        raise ValueError(
            f"asked for {start}-{end}, but the file only covers "
            f"{first.date()} to {last.date()} of that range"
        )
    return selected


def add_masks(
    data: pd.DataFrame,
    *,
    nee_column: str,
    sigma_column: str,
    qc_column: str,
    qc_threshold: float,
) -> pd.DataFrame:
    """Add True/False columns saying which days can enter the likelihood.

    No day is removed: the forward model must run through every day, and the
    likelihood only scores days where MASK_USE is True.
    """
    if not 0.0 <= qc_threshold <= 1.0:
        raise ValueError(f"qc_threshold must be between 0 and 1, got {qc_threshold}")

    out = data.copy()
    nee = out[nee_column]
    sigma = out[sigma_column]
    qc = out[qc_column]

    out[MASK_NEE] = nee.notna()
    out[MASK_QC] = qc.notna() & (qc >= qc_threshold)
    out[MASK_SIGMA] = sigma.notna() & (sigma > 0.0)
    out[MASK_USE] = out[MASK_NEE] & out[MASK_QC] & out[MASK_SIGMA]
    return out


def screening_report(data: pd.DataFrame) -> pd.DataFrame:
    """Sort every day into exactly one row: assimilable, or why it was excluded.

    The rows do not overlap, so they always sum to the number of days.
    """
    nee, qc, sigma = data[MASK_NEE], data[MASK_QC], data[MASK_SIGMA]

    groups = {
        "assimilable": nee & qc & sigma,
        "excluded: NEE missing": ~nee,
        "excluded: QC below threshold only": nee & ~qc & sigma,
        "excluded: RANDUNC missing or zero only": nee & qc & ~sigma,
        "excluded: both QC and RANDUNC": nee & ~qc & ~sigma,
    }

    total = len(data)
    report = pd.DataFrame({"days": {name: int(mask.sum()) for name, mask in groups.items()}})
    report["percent"] = (100.0 * report["days"] / total).round(1)

    if report["days"].sum() != total:
        raise AssertionError("report rows do not add up to the number of days")
    if report.loc["assimilable", "days"] != int(data[MASK_USE].sum()):
        raise AssertionError("report disagrees with the assimilable mask")
    return report


def load_screened_block(
    config: dict[str, Any] | None = None, *, block: str = "calibration"
) -> pd.DataFrame:
    """Load the file, keep one year block, and add the likelihood masks.

    Every day of the block is returned; use data[MASK_USE] for assimilable days.
    """
    if config is None:
        config = load_config()
    cols = config["columns"]

    data = load_fluxnet_dd(
        resolve_path(config["paths"]["fluxnet_dd"]),
        timestamp_column=cols["timestamp"],
        required_columns=[cols["nee"], cols["sigma"], cols["qc"]],
        missing_value=float(config["data"]["missing_value"]),
    )
    start, end = year_block(config, block)
    data = select_years(data, start, end)
    return add_masks(
        data,
        nee_column=cols["nee"],
        sigma_column=cols["sigma"],
        qc_column=cols["qc"],
        qc_threshold=float(config["data"]["qc_threshold"]),
    )


if __name__ == "__main__":
    cfg = load_config()
    print(f"QC threshold: {cfg['data']['qc_threshold']}")
    for name in ("calibration", "evaluation"):
        start, end = year_block(cfg, name)
        print(f"\n{name} {start}-{end}")
        print(screening_report(load_screened_block(cfg, block=name)).to_string())