"""Print the data columns the model uses, for a screenshot on the data slide.

Run from the repository root:
    python scripts/show_data_columns.py              # first days of 1997
    python scripts/show_data_columns.py --start 2005-06-28 --days 8

Shows the FLUXNET2015 columns used (NEE, its QC and RANDUNC, sunlight, CO2),
the daily highest and lowest temperature made from the half-hourly file (D18),
and whether each day enters the likelihood (D12). Nothing is saved.
"""

import argparse

import pandas as pd

from dalec.config import load_config
from dalec.data_io import MASK_USE, load_model_block


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--start", default=None, help="first date to show, YYYY-MM-DD")
    parser.add_argument("--days", type=int, default=8, help="number of days to show")
    args = parser.parse_args()

    config = load_config()
    data = load_model_block(config, block="calibration")
    cols = config["columns"]
    drivers = config["drivers"]["model"]
    extremes = config["drivers"]["temperature_extremes"]
    shown = {
        cols["nee"]: "NEE_VUT_REF",
        cols["qc"]: "NEE_VUT_REF_QC",
        cols["sigma"]: "NEE_VUT_REF_RANDUNC",
        drivers["sw_in"]: "SW_IN_F",
        drivers["co2"]: "CO2_F_MDS",
        extremes["t_max"]: "TA max (from HH)",
        extremes["t_min"]: "TA min (from HH)",
        MASK_USE: "used in fit",
    }
    table = data[list(shown)].rename(columns=shown)
    if args.start:
        table = table.loc[pd.Timestamp(args.start):]
    table = table.head(args.days)

    with pd.option_context("display.width", 140, "display.max_columns", 20,
                           "display.float_format", "{:.2f}".format):
        print(table.to_string())
    print(f"\n{len(data):,} days in 1997-2010, of which {int(data[MASK_USE].sum()):,} enter the fit "
          f"(QC at least {config['data']['qc_threshold']} and RANDUNC present)")


if __name__ == "__main__":
    main()
