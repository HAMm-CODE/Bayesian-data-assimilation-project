"""Test fitted models on the held-back years 2011-2014 and compare them (DECISIONS.md D61).

Nothing is sampled again. For every fit named on the command line, posterior
draws are run forward from 1 January 1997 straight through 31 December 2014, so
the pools at the start of 2011 come from the fitted years and the model has never
seen the test years. The rules below were fixed in D61 before the full fits were read.

Run from the repository root, for example:
    python scripts/06_evaluate.py fit_model_error_cold fit_model_error_dense_ta95

For each fit:
    0. the D49 computation checks; a fit that fails them is not tested or compared
    1. test years (decision): log predictive density on the 810 test days (QC at least
       0.75 and RANDUNC present, the fitting rule D12), coverage of the 90% predictive
       interval (overall and by season), RMSE and bias, NEE summed per test year;
       RMSE and bias also on the QC-only test days (secondary)
    2. how believable the fitted forest is over 1997-2010 (D40 Tier 2, peak LAI)
    3. PSIS-LOO on the 1997-2010 assimilable days (reported only; optimistic, because
       neighbouring days are correlated)
Then the comparison: score difference, its standard error widened for the
correlation between consecutive days, and the verdict under the D61 rule.

Saves (in figures/):
    evaluation_<fit>[_vs_<fit>...].csv      every measure for every fit
    fig_Y_<fit>[_vs_<fit>...]_test_years    monthly NEE in 2011-2014: observed and each model

Option --allow-unconverged runs fits that fail D49 anyway (for testing the code
only): every result is then marked NOT FOR REPORTING.
"""

import argparse
import warnings
from typing import Any

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import logsumexp

from dalec.config import load_config, resolve_path, year_block
from dalec.data_io import MASK_NEE, MASK_QC, MASK_USE, load_record_with_drivers, select_years
from dalec.model_numpy import drivers_from_frame, run
from dalec.plausibility import evaluate, run_and_summarise
from dalec.plotting import BLUE, INK, INK_SOFT, ORANGE, UNIT_FLUX, apply_style, heading, save_figure
from dalec.posterior import QUANTITIES, posterior_parameters, sampler_health
from dalec.priors import ACM_PRIORS, COLD_PRIORS

SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
MODEL_COLOURS = (BLUE, ORANGE, "#e87ba4", "#008300")   # same validated set as the chain colours


# ---------------------------------------------------------------------------
# Loading and the D49 checks
# ---------------------------------------------------------------------------

def load_fit(label: str, config: dict[str, Any]):
    path = resolve_path(config["paths"]["results"]) / f"{label}.nc"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; give the file name in results/ without .nc")
    return az.from_netcdf(path)


def sampled_names(post) -> list[str]:
    """Every sampled quantity (constants held fixed, D54, are left out)."""
    names = list(QUANTITIES) + list(COLD_PRIORS) + list(ACM_PRIORS) + ["sigma_m"]
    return [n for n in names if n in post.data_vars and float(np.ptp(np.asarray(post[n].values))) > 0]


def computation_checks(idata, config: dict[str, Any]) -> tuple[bool, dict[str, Any]]:
    """D49 criteria 1 to 3, recomputed from the saved fit."""
    diag = config["diagnostics"]
    summary = az.summary(idata, var_names=sampled_names(idata.posterior))
    health = sampler_health(idata)
    numbers = {
        "divergences": health["divergences"],
        "max_rhat": float(summary["r_hat"].max()),
        "min_ess": float(min(summary["ess_bulk"].min(), summary["ess_tail"].min())),
    }
    passed = (numbers["divergences"] == 0 and numbers["max_rhat"] < float(diag["max_rhat"])
              and numbers["min_ess"] >= float(diag["min_ess"]))
    return passed, numbers


# ---------------------------------------------------------------------------
# Running the posterior forward
# ---------------------------------------------------------------------------

def frost_threshold(label: str, config: dict[str, Any]) -> float | None:
    frost = config["model"]["frost_cutoff"]
    on = "_frost" in label or bool(frost["enabled"])
    return float(frost["threshold_degc"]) if on else None


def forward(idata, data: pd.DataFrame, config: dict[str, Any], label: str, n: int, seed: int):
    """Daily model NEE (draws x days) and sigma_m per draw, for the same draws every time."""
    params, sigma_m, cold, a7 = posterior_parameters(idata, n, seed)
    drivers = drivers_from_frame(data, config)
    latitude = float(config["site"]["latitude_deg"])
    batch_size = int(config["prior_predictive"]["batch_size"])
    arrays = params.as_arrays()
    nee = []
    for start in range(0, sigma_m.size, batch_size):
        stop = min(start + batch_size, sigma_m.size)
        batch = type(params)(**{k: v[start:stop] for k, v in arrays.items()})
        with np.errstate(all="ignore"):
            out = run(batch, drivers, latitude_deg=latitude, frost_threshold_degc=frost_threshold(label, config),
                      cold=None if cold is None else (cold[0][start:stop], cold[1][start:stop]),
                      a7=None if a7 is None else a7[start:stop])
        nee.append(out.nee)
    return np.concatenate(nee), sigma_m


def believability(idata, data_cal: pd.DataFrame, config: dict[str, Any], label: str,
                  n: int, seed: int) -> dict[str, float]:
    """D40 Tier 2 share and peak LAI over 1997-2010, for the same draws."""
    params, _, cold, a7 = posterior_parameters(idata, n, seed)
    summary = run_and_summarise(params, drivers_from_frame(data_cal, config),
                                float(config["site"]["latitude_deg"]), frost_threshold(label, config),
                                data_cal.index.year.to_numpy(), None,
                                int(config["prior_predictive"]["batch_size"]),
                                verbose=False, cold=cold, a7=a7)
    checks = evaluate(summary)
    return {"tier2_share": float(checks["tier2_all"].mean()),
            "peak_lai_median": float(np.median(summary["peak_lai"])),
            "gpp_median": float(np.median(summary["gpp"])),
            "rh_median": float(np.median(summary["rh"]))}


# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------

def pointwise_loglik(model_nee: np.ndarray, observed: np.ndarray, sigma: np.ndarray,
                     sigma_m: np.ndarray) -> np.ndarray:
    """log Normal(observed | model, sqrt(RANDUNC^2 + sigma_m^2)) per draw and day (L1)."""
    sd = np.sqrt(sigma[None, :] ** 2 + sigma_m[:, None] ** 2)
    return stats.norm.logpdf(observed[None, :], model_nee, sd)


def lpd_per_day(loglik: np.ndarray) -> np.ndarray:
    """Log of the average predictive density over draws, one value per day."""
    return logsumexp(loglik, axis=0) - np.log(loglik.shape[0])


def test_measures(model_nee, sigma_m, observed, sigma, index, use, qc_only, interval, seed):
    """Every D61 measure on the test days."""
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    lpd = lpd_per_day(pointwise_loglik(model_nee[:, use], observed[use], sigma[use], sigma_m))
    sd = np.sqrt(sigma[use][None, :] ** 2 + sigma_m[:, None] ** 2)
    noise = np.random.default_rng(seed).standard_normal(sd.shape)
    replicated = model_nee[:, use] + sd * noise
    low, high = np.percentile(replicated, [lo, hi], axis=0)
    inside = (observed[use] >= low) & (observed[use] <= high)
    mean_model = model_nee.mean(axis=0)
    months = index.month.to_numpy()[use]

    out: dict[str, Any] = {
        "test days": int(use.sum()),
        "log predictive density": float(lpd.sum()),
        "coverage, all": float(inside.mean()),
    }
    for season, members in SEASONS.items():
        out[f"coverage, {season}"] = float(inside[np.isin(months, members)].mean())
    for name, mask in (("", use), (" (QC only)", qc_only)):
        residual = observed[mask] - mean_model[mask]
        out[f"RMSE{name}"] = float(np.sqrt(np.mean(residual ** 2)))
        out[f"bias, observed minus model{name}"] = float(residual.mean())
        out[f"days{name}"] = int(mask.sum())
    years = index.year.to_numpy()
    for y in np.unique(years[use]):
        days = use & (years == y)
        out[f"{y} NEE sum, observed"] = float(observed[days].sum())
        out[f"{y} NEE sum, model median"] = float(np.median(model_nee[:, days].sum(axis=1)))
    return out, lpd


def psis_loo(loglik: np.ndarray) -> tuple[dict[str, float], np.ndarray]:
    """PSIS-LOO (Vehtari et al. 2017) from a draws x days log-likelihood matrix."""
    idata = az.from_dict({"posterior": {"draw_index": np.arange(loglik.shape[0], dtype=float)[None, :]},
                          "log_likelihood": {"nee": loglik[None, :, :]}})
    with warnings.catch_warnings():          # high Pareto k is counted and reported below
        warnings.simplefilter("ignore", UserWarning)
        result = az.loo(idata, pointwise=True)
    k = np.asarray(result.pareto_k).reshape(-1)
    return ({"PSIS-LOO elpd (1997-2010, optimistic)": float(result.elpd),
             "PSIS-LOO standard error": float(result.se),
             "Pareto k above the good limit": int((k > float(result.good_k)).sum()),
             "Pareto k above 1": int((k > 1).sum())},
            np.asarray(result.elpd_i).reshape(-1))


def widened_se(diff: np.ndarray, dates: pd.DatetimeIndex) -> tuple[float, float, float]:
    """Standard error of a summed day-by-day difference: plain (Vehtari et al. 2017) and
    widened by sqrt((1 + r) / (1 - r)) for the correlation r between consecutive days (D61)."""
    plain = float(np.sqrt(diff.size) * diff.std(ddof=1))
    consecutive = np.diff(dates.values).astype("timedelta64[D]").astype(int) == 1
    r = float(np.corrcoef(diff[:-1][consecutive], diff[1:][consecutive])[0, 1]) if consecutive.sum() > 2 else 0.0
    factor = np.sqrt((1 + r) / (1 - r)) if 0 < r < 1 else 1.0
    return plain, r, plain * factor


# ---------------------------------------------------------------------------
# Figure Y
# ---------------------------------------------------------------------------

def test_years_figure(results, observed, index, use, interval, config, name) -> None:
    """Monthly mean NEE over the test days: observed (dots) and each model (median, band)."""
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 12, "xtick.labelsize": 10.5,
                         "ytick.labelsize": 10.5, "legend.fontsize": 11})
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    month_id = index.to_period("M")
    months = month_id[use].unique()
    centre = np.array([m.to_timestamp() + pd.Timedelta(days=14) for m in months])
    fig, ax = plt.subplots(figsize=(12, 4.6))
    for k, (label, res) in enumerate(results.items()):
        colour = MODEL_COLOURS[k % len(MODEL_COLOURS)]
        monthly = np.stack([res["nee"][:, use & (month_id == m)].mean(axis=1) for m in months], axis=1)
        low, mid, high = np.percentile(monthly, [lo, 50, hi], axis=0)
        ax.fill_between(centre, low, high, color=colour, alpha=0.18, lw=0)
        ax.plot(centre, mid, color=colour, lw=2, label=f"{label}: median and {interval:.0%} band")
    obs = np.array([observed[use & (month_id == m)].mean() for m in months])
    ax.plot(centre, obs, "o", color=INK, ms=4, label="observed (test days)")
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    ax.set_ylabel(f"monthly mean NEE\n({UNIT_FLUX})")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncols=len(results) + 1, frameon=False)
    fig.tight_layout()
    save_figure(fig, name, config)


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("labels", nargs="+", help="fits in results/, without .nc")
    parser.add_argument("--allow-unconverged", action="store_true",
                        help="testing only: also run fits that fail D49 (results NOT FOR REPORTING)")
    args = parser.parse_args()

    config = load_config()
    apply_style()
    ev = config["evaluation"]
    n, seed, interval = int(ev["n_draws"]), int(ev["seed"]), float(ev["interval"])
    cal_start, cal_end = year_block(config, "calibration")
    test_start, test_end = year_block(config, "evaluation")
    data = select_years(load_record_with_drivers(config), cal_start, test_end)
    data_cal = select_years(data, cal_start, cal_end)
    cols = config["columns"]
    observed = data[cols["nee"]].to_numpy(dtype=float)
    sigma = data[cols["sigma"]].to_numpy(dtype=float)
    years = data.index.year.to_numpy()
    in_test = (years >= test_start) & (years <= test_end)
    in_cal = (years >= cal_start) & (years <= cal_end)
    use = data[MASK_USE].to_numpy(dtype=bool)
    test_use = in_test & use
    test_qc_only = in_test & data[MASK_NEE].to_numpy(dtype=bool) & data[MASK_QC].to_numpy(dtype=bool)
    cal_use = in_cal & use

    results: dict[str, dict[str, Any]] = {}
    rows = {}
    for label in args.labels:
        heading(f"{label}")
        idata = load_fit(label, config)
        passed, numbers = computation_checks(idata, config)
        print(f"  D49: divergences {numbers['divergences']}, largest R-hat {numbers['max_rhat']:.4f}, "
              f"smallest ESS {numbers['min_ess']:.0f}  ->  {'PASS' if passed else 'FAIL'}")
        if not passed and not args.allow_unconverged:
            print("  not tested or compared (D61: only fits that pass D49)")
            continue
        model_nee, sigma_m = forward(idata, data, config, label, n, seed)
        measures, lpd = test_measures(model_nee, sigma_m, observed, sigma, data.index,
                                      test_use, test_qc_only, interval, seed + 1)
        fitted = cal_use.copy()             # PSIS-LOO only on the days this fit used
        if label.endswith("_quick"):
            fitted &= years < cal_start + int(config["sampler"]["quick"]["years"])
        loo, elpd_i = psis_loo(pointwise_loglik(model_nee[:, fitted], observed[fitted],
                                                sigma[fitted], sigma_m))
        forest = believability(idata, data_cal, config, label, n, seed)
        row = {"passes D49": passed, **numbers, **measures, **loo,
               "Tier 2 share (1997-2010)": forest["tier2_share"],
               "peak LAI, median (site 2.5 to 4)": forest["peak_lai_median"],
               "annual GPP, median": forest["gpp_median"], "annual Rh, median": forest["rh_median"]}
        if not passed:
            row["NOTE"] = "NOT FOR REPORTING: fails D49"
        rows[label] = row
        results[label] = {"nee": model_nee, "lpd": lpd, "elpd_i": elpd_i, "row": row}
        with pd.option_context("display.float_format", "{:.4g}".format):
            print(pd.Series(row).to_string())

    if not results:
        print("\nNo fit passed D49: nothing to test or compare.")
        return

    table = pd.DataFrame(rows)
    if len(results) > 1:
        heading("Comparison on the test years (D61)")
        ranked = sorted(results, key=lambda k: results[k]["row"]["log predictive density"], reverse=True)
        best = ranked[0]
        margin = float(ev["clear_margin"])
        dates = data.index[test_use]
        for other in ranked[1:]:
            diff = results[best]["lpd"] - results[other]["lpd"]
            plain, r, wide = widened_se(diff, dates)
            clear = diff.sum() > margin * wide
            loo_diff = results[best]["elpd_i"] - results[other]["elpd_i"]
            print(f"  {best} minus {other}:")
            print(f"    test-years score difference {diff.sum():.1f}; standard error {plain:.1f} "
                  f"(widened {wide:.1f}, correlation of consecutive days {r:+.2f})")
            print(f"    verdict: {'CLEARLY BETTER: ' + best if clear else 'NO CLEAR DIFFERENCE'} "
                  f"(rule: difference above {margin:g} x widened standard error = {margin * wide:.1f})")
            if not clear:
                print("    -> D61: report both; prefer the model whose forest is believable "
                      "(Tier 2 share and peak LAI in the table)")
            print(f"    PSIS-LOO 1997-2010 (reported only, optimistic): difference {loo_diff.sum():.1f}, "
                  f"standard error {np.sqrt(loo_diff.size) * loo_diff.std(ddof=1):.1f}")
            table.loc[f"vs {other}: score difference", best] = float(diff.sum())
            table.loc[f"vs {other}: widened standard error", best] = wide
            table.loc[f"vs {other}: verdict", best] = "clearly better" if clear else "no clear difference"

    name = "_vs_".join(args.labels if len(args.labels) <= 3 else args.labels[:3])
    out_dir = resolve_path(config["paths"]["figures"])
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / f"evaluation_{name}.csv")
    print(f"\n  saved {out_dir / f'evaluation_{name}.csv'}")
    test_years_figure(results, observed, data.index, test_use, interval, config,
                      f"fig_Y_{name}_test_years")
    if any(not res["row"]["passes D49"] for res in results.values()):
        print("\n  NOT FOR REPORTING: at least one fit fails D49 (--allow-unconverged)")


if __name__ == "__main__":
    main()
