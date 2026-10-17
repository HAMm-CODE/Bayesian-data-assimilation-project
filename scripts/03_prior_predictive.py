"""Prior predictive check (Gelman et al. 2020, Section 2.4; DECISIONS.md D40).

Draws parameter sets from all priors at once (D34 to D39), runs DALEC2 with the
real 1997-2010 weather for each, and asks: do the priors describe a believable
boreal pine forest? The pass criteria were fixed in D40 before this was run.
No NEE observations are used in the check; the observed NEE cycle appears in
Figure L for context only. The D40 checks themselves live in
`src/dalec/plausibility.py`, shared with the fake-data check (script 04).

Run from the repository root:
    python scripts/03_prior_predictive.py          main model (DALEC2)
    python scripts/03_prior_predictive.py --cold   cold-limit model (D55) with its priors (D57)

Makes two figures (names end in _cold for the cold-limit model):
    L  fig_L_prior_predictive_nee      NEE through the year: prior range vs observed (context)
    M  fig_M_prior_predictive_checks   the D40 quantities, with their plausible bands
"""

import argparse

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path, year_block
from dalec.data_io import MASK_NEE, MASK_QC, load_model_block
from dalec.model_numpy import drivers_from_frame
from dalec.plausibility import (
    TIER1_TARGET, TIER2_BANDS, TIER2_TARGET, evaluate, run_and_summarise,
)
from dalec.plotting import (
    BLUE, GREY_LIGHT, INK, INK_SOFT, UNIT_FLUX,
    apply_style, climatology, day_of_year_365, heading, save_figure, set_year_axis,
)
from dalec.priors import (
    COLD_MODEL_PRIOR_CHANGES, COLD_PRIORS, PRIORS, edc_log_weight, priors_for, sample_cold,
    sample_prior, to_prior_scale,
)

# ---------------------------------------------------------------------------
# Running the draws
# ---------------------------------------------------------------------------

def doy_weights(index: pd.DatetimeIndex) -> np.ndarray:
    """Matrix that turns a daily series into its mean for each day of the year (1..365)."""
    doy = day_of_year_365(index).to_numpy()
    weights = np.zeros((len(index), 365))
    keep = ~np.isnan(doy)
    weights[np.flatnonzero(keep), doy[keep].astype(int) - 1] = 1.0
    return weights / weights.sum(axis=0)


# ---------------------------------------------------------------------------
# Applying D40
# ---------------------------------------------------------------------------

def report(result: dict[str, Any], title: str) -> None:
    heading(title)
    n = result["tier1_fail"].size
    t1 = result["tier1_fail"].mean()
    print(f"  Tier 1, physically impossible: {t1:6.1%} of {n} draws "
          f"(target at most {TIER1_TARGET:.0%})  -> {'PASS' if t1 <= TIER1_TARGET else 'FAIL'}")
    for reason, flag in result["tier1_reasons"].items():
        if flag.any():
            print(f"      {reason:24s} {flag.mean():6.1%}")
    t2 = result["tier2_all"].mean()
    print(f"  Tier 2, plausible for Hyytiala: {t2:6.1%} pass every check "
          f"(target at least {TIER2_TARGET:.0%})  -> {'PASS' if t2 >= TIER2_TARGET else 'FAIL'}")
    for name, flag in result["tier2"].items():
        print(f"      {name:14s} {flag.mean():6.1%} inside the band")


def drivers_of_failure(raw: dict[str, np.ndarray], failed: np.ndarray, label: str) -> dict[str, float]:
    """Which priors push draws into failure? For every sampled quantity, the
    standardised difference of means between failing and passing draws, on the
    prior's own scale. 0 = no link; about 0.2 or more (either sign) = worth a look."""
    if failed.all() or not failed.any():
        return {}
    links = {}
    for name, values in raw.items():
        if name == "theta_min":      # follows from theta_lit and h
            continue
        x = to_prior_scale(name, values)
        links[name] = float((x[failed].mean() - x[~failed].mean()) / x.std())
    print(f"\n  All sampled quantities, linked to '{label}' (failing minus passing):")
    for name, diff in sorted(links.items(), key=lambda kv: -abs(kv[1])):
        direction = "higher" if diff > 0 else "lower"
        flag = "  <-- worth a look" if abs(diff) >= 0.2 else ""
        print(f"      {name:22s} {diff:+.2f}  (failing draws have {direction} values){flag}")
    return links


def parameter_table(raw, tier1_links, tier2_links, config, suffix: str = "") -> None:
    """One row per sampled quantity: prior spread and link to each failure. Saved as csv."""
    rows = []
    for name, values in raw.items():
        rows.append({
            "parameter": name,
            "p05": np.percentile(values, 5),
            "median": np.median(values),
            "p95": np.percentile(values, 95),
            "link_tier1_failure": tier1_links.get(name, np.nan),
            "link_tier2_failure": tier2_links.get(name, np.nan),
        })
    table = pd.DataFrame(rows).set_index("parameter")
    path = resolve_path(config["paths"]["figures"]) / f"prior_predictive_parameters{suffix}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, float_format="%.6g")
    print(f"\n  saved {path}")


def edc_weights(m: dict[str, np.ndarray]) -> np.ndarray:
    """Weight of each draw under the soft steady-state constraint on the foliage (D44), max = 1."""
    log_w = edc_log_weight(m["foliage_rs"])
    return np.exp(log_w - np.max(log_w[np.isfinite(log_w)]))


def weighted_report(result: dict[str, Any], w: np.ndarray, title: str) -> None:
    """D40 rates when every draw counts with its soft-constraint weight (D44).

    This is importance re-weighting (Gelman et al. 2020, Section 6.3): it shows
    what the prior with the soft constraint would give, without new draws.
    """
    heading(title)
    total = w.sum()
    ess = total**2 / np.sum(w**2)
    print(f"  effective number of draws after weighting: {ess:.0f} of {w.size}")
    print(f"  draws with almost no weight (below 1%): {np.mean(w < 0.01):.1%}")
    t1 = np.sum(w * result["tier1_fail"]) / total
    print(f"  Tier 1, physically impossible: {t1:6.1%} (target at most {TIER1_TARGET:.0%})"
          f"  -> {'PASS' if t1 <= TIER1_TARGET else 'FAIL'}")
    for reason, flag in result["tier1_reasons"].items():
        if flag.any():
            print(f"      {reason:24s} {np.sum(w * flag) / total:6.1%}")
    t2 = np.sum(w * result["tier2_all"]) / total
    print(f"  Tier 2, plausible for Hyytiala: {t2:6.1%} pass every check (target at least "
          f"{TIER2_TARGET:.0%})  -> {'PASS' if t2 >= TIER2_TARGET else 'FAIL'}")
    for name, flag in result["tier2"].items():
        print(f"      {name:14s} {np.sum(w * flag) / total:6.1%} inside the band")
    return t1, t2


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def figure_l(m: dict[str, np.ndarray], observed: pd.DataFrame, config: dict[str, Any],
             suffix: str = "") -> None:
    """Prior predictive range of the NEE seasonal cycle, observed cycle for context."""
    ok = ~evaluate(m)["tier1_fail"]
    cycle = m["nee_cycle"][ok]
    days = np.arange(1, 366)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.fill_between(days, np.percentile(cycle, 5, axis=0), np.percentile(cycle, 95, axis=0),
                    color=BLUE, alpha=0.2, lw=0, label="prior predictive, 5th to 95th percentile")
    ax.fill_between(days, np.percentile(cycle, 25, axis=0), np.percentile(cycle, 75, axis=0),
                    color=BLUE, alpha=0.35, lw=0, label="prior predictive, middle half")
    ax.plot(days, np.median(cycle, axis=0), color=BLUE, lw=1.5, label="prior predictive median")
    ax.plot(observed.index, observed["median"], color=INK, lw=1.2, ls="--",
            label="observed median (context only, not used in the check)")
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    set_year_axis(ax)
    start, end = year_block(config, "calibration")
    model = " (cold-limit model)" if suffix == "_cold" else ""
    ax.set_title(f"NEE through the year implied by the priors{model}, {start}-{end}")
    ax.set_ylabel(f"NEE ({UNIT_FLUX})\nnegative = forest takes up carbon")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncols=2)
    save_figure(fig, f"fig_L_prior_predictive_nee{suffix}", config)


def figure_m(m: dict[str, np.ndarray], config: dict[str, Any], suffix: str = "") -> None:
    """Histograms of the D40 quantities with their plausible bands shaded."""
    ok = ~evaluate(m)["tier1_fail"]
    fig, axes = plt.subplots(2, 4, figsize=(13, 6))
    for ax, (name, (low, high, label)) in zip(axes.flat, TIER2_BANDS.items()):
        values = m[name][ok]
        lo, hi = np.percentile(values, [0.5, 99.5])
        span_lo = lo if np.isinf(low) else min(lo, low)
        span_hi = hi if np.isinf(high) else max(hi, high)
        ax.axvspan(low if np.isfinite(low) else span_lo, high if np.isfinite(high) else span_hi,
                   color=GREY_LIGHT, lw=0, zorder=0)
        ax.hist(np.clip(values, span_lo, span_hi), bins=40, color=BLUE, alpha=0.8, zorder=2)
        ax.set_xlabel(label, fontsize=8)
        ax.set_yticks([])
        ax.grid(False)
    ax = axes.flat[-1]
    ax.axis("off")
    ax.text(0, 0.6, "grey band = plausible range (D40)\nbars = draws that are\nphysically possible",
            fontsize=9, va="center")
    model = ", cold-limit model" if suffix == "_cold" else ""
    fig.suptitle(f"What the priors imply about the forest (prior predictive, D40{model})", fontsize=10)
    fig.tight_layout()
    save_figure(fig, f"fig_M_prior_predictive_checks{suffix}", config)


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cold", action="store_true",
                        help="cold-limit model (D55) with its own priors (D57)")
    args = parser.parse_args()

    config = load_config()
    apply_style()
    settings = config["prior_predictive"]

    data = load_model_block(config, block="calibration")
    drivers = drivers_from_frame(data, config)
    years = data.index.year.to_numpy()
    weights = doy_weights(pd.DatetimeIndex(data.index))
    latitude = float(config["site"]["latitude_deg"])
    batch_size = int(settings["batch_size"])
    qc_pass = data[MASK_NEE] & data[MASK_QC]
    observed = climatology(data[config["columns"]["nee"]].where(qc_pass))

    rng = np.random.default_rng(int(settings["seed"]))
    if args.cold:
        cold_model(drivers, years, weights, latitude, batch_size, observed, settings, config, rng)
        return
    params, raw = sample_prior(int(settings["n_draws"]), rng)

    heading(f"Running {settings['n_draws']} prior draws over {drivers.n_days} days "
            "(main model: no frost cut-off)")
    main_run = run_and_summarise(params, drivers, latitude, None, years, weights, batch_size)
    main_result = evaluate(main_run)
    report(main_result, "Prior predictive check, main model (D40)")
    tier1_links = drivers_of_failure(raw, main_result["tier1_fail"], "Tier 1 failure")
    tier2_links = drivers_of_failure(raw, ~main_result["tier2_all"], "Tier 2 failure")
    parameter_table(raw, tier1_links, tier2_links, config)

    heading("Same draws with the frost switch on (D24), reported alongside")
    frost = float(config["model"]["frost_cutoff"]["threshold_degc"])
    frost_run = run_and_summarise(params, drivers, latitude, frost, years, weights, batch_size)
    report(evaluate(frost_run), f"Prior predictive check, frost cut-off at {frost} degC")

    weighted_report(main_result, edc_weights(main_run),
                    "With the soft foliage constraint (D44), main model, by re-weighting")
    weighted_report(evaluate(frost_run), edc_weights(frost_run),
                    "With the soft foliage constraint (D44), frost switch on, by re-weighting")

    heading("Medians of the plausible draws (main model)")
    ok = ~main_result["tier1_fail"]
    for name in TIER2_BANDS:
        print(f"  {name:14s} median {np.median(main_run[name][ok]):9.2f}")

    figure_l(main_run, observed, config)
    figure_m(main_run, config)


# The main model was accepted with Tier 1 at 6.4% under the soft foliage
# constraint (F15, a narrow miss of the 5% target). The cold-limit model may not
# do clearly worse: at most 7% (D57; with about 800 effective draws the Monte
# Carlo error of such a rate is about 0.9 percentage points).
COLD_MODEL_TIER1_LIMIT = 0.07


def cold_model(drivers, years, weights, latitude, batch_size, observed, settings, config, rng) -> None:
    """Prior predictive check of the cold-limit model (D55) with its priors (D57).

    The same draws are also run with the `ceff` prior of the main model (D42), by
    moving each `ceff` draw to the same position within the D42 prior, so the two
    differ only in `ceff`.
    """
    n = int(settings["n_draws"])
    priors = priors_for(cold=True)
    params, raw = sample_prior(n, rng, priors)
    raw.update(sample_cold(n, rng))
    thresholds = tuple(raw[name] for name in COLD_PRIORS)
    changed = ", ".join(f"{k} {v}" for k, v in COLD_MODEL_PRIOR_CHANGES.items())

    heading(f"Running {n} prior draws over {drivers.n_days} days "
            f"(cold-limit model, D55; priors changed for it, D57: {changed})")
    cold_run = run_and_summarise(params, drivers, latitude, None, years, weights, batch_size,
                                 cold=thresholds)
    cold_result = evaluate(cold_run)
    report(cold_result, "Prior predictive check, cold-limit model (D40, D57)")
    tier1_links = drivers_of_failure(raw, cold_result["tier1_fail"], "Tier 1 failure")
    tier2_links = drivers_of_failure(raw, ~cold_result["tier2_all"], "Tier 2 failure")
    parameter_table(raw, tier1_links, tier2_links, config, suffix="_cold")
    t1, t2 = weighted_report(cold_result, edc_weights(cold_run),
                             "With the soft foliage constraint (D44), cold-limit model, by re-weighting")

    # Same draws, ceff moved to the same position within the main model's prior (D42).
    _, median_new, sd_new = priors["ceff"]
    _, median_old, sd_old = PRIORS["ceff"]
    z = np.log(raw["ceff"] / median_new) / sd_new
    arrays = params.as_arrays()
    arrays["ceff"] = median_old * np.exp(sd_old * z)
    old_run = run_and_summarise(type(params)(**arrays), drivers, latitude, None, years, weights,
                                batch_size, verbose=False, cold=thresholds)
    old_result = evaluate(old_run)
    report(old_result, "For comparison: same draws, ceff prior of the main model (D42)")
    t1_old, t2_old = weighted_report(old_result, edc_weights(old_run),
                                     "For comparison, with the soft foliage constraint (D44)")

    heading("Medians of the plausible draws (cold-limit model; D42 ceff in brackets)")
    ok, ok_old = ~cold_result["tier1_fail"], ~old_result["tier1_fail"]
    for name in TIER2_BANDS:
        print(f"  {name:14s} median {np.median(cold_run[name][ok]):9.2f}"
              f"   ({np.median(old_run[name][ok_old]):9.2f})")

    heading("D57 rule, fixed before running (with the soft foliage constraint, as fitted)")
    rule = {
        f"Tier 2 at least {TIER2_TARGET:.0%}": t2 >= TIER2_TARGET,
        f"Tier 1 at most {COLD_MODEL_TIER1_LIMIT:.0%} (main model accepted at 6.4%, F15)":
            t1 <= COLD_MODEL_TIER1_LIMIT,
    }
    for text, ok_rule in rule.items():
        print(f"  {'PASS' if ok_rule else 'FAIL'}  {text}")
    print(f"  (cold model: Tier 1 {t1:.1%}, Tier 2 {t2:.1%}; "
          f"with the D42 ceff prior: Tier 1 {t1_old:.1%}, Tier 2 {t2_old:.1%})")
    verdict = all(rule.values())
    print(f"\n  -> {'use these priors for the cold-limit model' if verdict else 'keep the D42 ceff prior for the cold-limit model (D57 option B)'}")

    figure_l(cold_run, observed, config, suffix="_cold")
    figure_m(cold_run, config, suffix="_cold")


if __name__ == "__main__":
    main()