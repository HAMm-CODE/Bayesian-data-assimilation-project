"""What the fake data of a fake-data check look like (D48), as a figure and a csv.

Rebuilds, from a saved fake-data check, the three series it rests on:
    true model NEE  DALEC2 run with the true forest (the truth column of the recovery csv)
    fake NEE        the true model NEE plus noise of the real size (RANDUNC), on the real
                    assimilable days only; read back from the saved run, exactly as fitted
    real NEE        the measured NEE on the same days, for comparison
Nothing is sampled.

Run from the repository root with the run name (results/<name>.nc):
    python scripts/09_fake_data_figure.py fake_data_check_cold_dense_ta99 [--year 2003]

Saves:
    figures/fake_data_series_<name>.csv    one row per day: date, used, real, RANDUNC, true model, fake
    Z  fig_Z_fake_data_<name>              top: monthly means 1997-2010; bottom: one year, daily
"""

import argparse

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.model_numpy import Parameters, drivers_from_frame, run
from dalec.plotting import BLUE, GREY, INK, ORANGE, UNIT_FLUX, apply_style, save_figure
from dalec.priors import COLD_PRIORS, FRACTIONS, shares_from_fractions


def true_forest(table: pd.DataFrame):
    """Parameters and cold thresholds of the true forest, from the recovery table."""
    truth = table["truth"].astype(float).to_dict()
    truth.update(shares_from_fractions(*(truth[k] for k in ("f_auto", *FRACTIONS))))
    params = Parameters(**{k: truth[k] for k in Parameters.__dataclass_fields__})
    cold = (tuple(truth[k] for k in COLD_PRIORS) if all(k in truth for k in COLD_PRIORS) else None)
    return params, cold


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run", help="fake-data check name, without .nc")
    parser.add_argument("--year", type=int, default=2003, help="the year shown day by day")
    args = parser.parse_args()
    config = load_config()
    apply_style()

    suffix = args.run.removeprefix("fake_data_check")
    table = pd.read_csv(resolve_path(config["paths"]["figures"]) / f"fake_data_recovery{suffix}.csv", index_col=0)
    idata = az.from_netcdf(resolve_path(config["paths"]["results"]) / f"{args.run}.nc")
    fake_values = np.asarray(idata.observed_data["nee"].values, dtype=float)

    data = load_model_block(config, block="calibration")
    use = data[MASK_USE].to_numpy(dtype=bool)
    n_days = int(np.searchsorted(np.cumsum(use), fake_values.size) + 1)   # quick runs cover fewer days
    data, use = data.iloc[:n_days], use[:n_days]
    drivers = drivers_from_frame(data, config)
    params, cold = true_forest(table)
    true_nee = run(params, drivers, latitude_deg=float(config["site"]["latitude_deg"]), cold=cold).nee
    fake = np.full(n_days, np.nan)
    fake[use] = fake_values
    real = data[config["columns"]["nee"]].to_numpy(dtype=float)
    sigma = data[config["columns"]["sigma"]].to_numpy(dtype=float)

    series = pd.DataFrame({"used": use, "real_nee": np.where(use, real, np.nan), "randunc": np.where(use, sigma, np.nan),
                           "true_model_nee": true_nee, "fake_nee": fake}, index=data.index)
    series.index.name = "date"
    out = resolve_path(config["paths"]["figures"]) / f"fake_data_series{suffix}.csv"
    series.to_csv(out, float_format="%.4f")
    print(f"  saved {out}")
    check = (series["fake_nee"] - series["true_model_nee"])[use] / series["randunc"][use]
    print(f"  {int(use.sum())} fake days; (fake - true model) / RANDUNC: mean {check.mean():+.3f}, sd {check.std():.3f} "
          "(about 0 and 1 if the noise is as designed)")

    fig, (top, bottom) = plt.subplots(2, 1, figsize=(12, 7.5))
    monthly = series.resample("MS").mean()
    top.plot(monthly.index, monthly["true_model_nee"], color=INK, lw=2, label="true model NEE (all days)")
    top.plot(monthly.index, monthly["fake_nee"], color=ORANGE, lw=1.6, label="fake NEE (assimilable days)")
    top.plot(monthly.index, monthly["real_nee"], color=GREY, lw=1.6, label="measured NEE (same days)")
    top.axhline(0, color=GREY, lw=0.8)
    top.set_ylabel(f"monthly mean NEE ({UNIT_FLUX})")
    top.legend(loc="lower left", frameon=False, ncol=3)

    year = series[series.index.year == args.year]
    if year.empty:
        year = series[series.index.year == series.index.year[0]]
    bottom.plot(year.index, year["true_model_nee"], color=INK, lw=2, label="true model NEE")
    bottom.errorbar(year.index, year["fake_nee"], yerr=year["randunc"], fmt="o", ms=3.5, color=ORANGE,
                    ecolor=ORANGE, elinewidth=0.6, alpha=0.8, label="fake NEE ± RANDUNC")
    bottom.plot(year.index, year["real_nee"], "o", ms=3, color=BLUE, alpha=0.5, label="measured NEE")
    bottom.axhline(0, color=GREY, lw=0.8)
    bottom.set_ylabel(f"daily NEE ({UNIT_FLUX})")
    bottom.set_title(f"{year.index.year[0]}, day by day", loc="left")
    bottom.legend(loc="lower left", frameon=False, ncol=3)
    fig.tight_layout()
    save_figure(fig, f"fig_Z_fake_data{suffix}", config)


if __name__ == "__main__":
    main()
