"""Prior predictive check (Gelman et al. 2020, Section 2.4; DECISIONS.md D40).

Draws parameter sets from all priors at once (D34 to D39), runs DALEC2 with the
real 1997-2010 weather for each, and asks: do the priors describe a believable
boreal pine forest? The pass criteria were fixed in D40 before this was run.
No NEE observations are used in the check; the observed NEE cycle appears in
Figure L for context only.

Run from the repository root:
    python scripts/03_prior_predictive.py

Makes two figures:
    L  fig_L_prior_predictive_nee      NEE through the year: prior range vs observed (context)
    M  fig_M_prior_predictive_checks   the D40 quantities, with their plausible bands
"""

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path, year_block
from dalec.data_io import MASK_NEE, MASK_QC, load_model_block
from dalec.model_numpy import POOLS, drivers_from_frame, run
from dalec.plotting import (
    BLUE, GREY_LIGHT, INK, INK_SOFT, UNIT_FLUX,
    apply_style, climatology, day_of_year_365, heading, save_figure, set_year_axis,
)
from dalec.priors import PRIORS, edc_log_weight, sample_prior, steady_state_ratio

# ---------------------------------------------------------------------------
# D40 criteria (copied from DECISIONS.md; change only together with D40)
# ---------------------------------------------------------------------------

TIER1_TARGET = 0.05   # at most 5% of draws physically impossible
TIER2_TARGET = 0.50   # at least half of the draws pass every Tier-2 check

TIER2_BANDS = {       # name: (low, high, label)
    "gpp": (500.0, 2200.0, "annual GPP (g C m$^{-2}$ yr$^{-1}$)"),
    "rh": (145.0, 740.0, "annual heterotrophic respiration (g C m$^{-2}$ yr$^{-1}$)"),
    "nee": (-700.0, 300.0, "annual NEE (g C m$^{-2}$ yr$^{-1}$)"),
    "peak_lai": (1.5, 8.0, "peak LAI (mean of yearly maxima)"),
    "foliage_ratio": (1.0 / 3.0, 3.0, "foliage pool, end / start"),
    "wood_ratio": (0.8, np.inf, "wood pool, end / start"),
    "soil_ratio": (0.7, 1.3, "soil pool, end / start"),
}


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


def summarise_batch(out, years: np.ndarray, weights: np.ndarray) -> dict[str, np.ndarray]:
    """Per-draw numbers needed by D40, from one batch of model runs."""
    year_values = np.unique(years)
    def annual(series):   # (draws, years) of yearly totals
        return np.stack([series[:, years == y].sum(axis=1) for y in year_values], axis=1)

    gpp_year = annual(out.gpp)
    lai_year_max = np.stack([out.lai[:, years == y].max(axis=1) for y in year_values], axis=1)
    finite = np.isfinite(out.nee).all(axis=1) & np.all(
        [np.isfinite(out.pools_end[p]) for p in POOLS], axis=0)
    return {
        "finite": finite,
        "negative_pool": np.any([out.pools_min[p] < 0 for p in POOLS], axis=0),
        "lai_max": out.lai.max(axis=1),
        "gpp_daily_max": out.gpp.max(axis=1),
        "gpp_annual_max": gpp_year.max(axis=1),
        "gpp": gpp_year.mean(axis=1),
        "rh": annual(out.rh).mean(axis=1),
        "nee": annual(out.nee).mean(axis=1),
        "peak_lai": lai_year_max.mean(axis=1),
        "foliage_ratio": out.pools_end["c_fol"] / out.pools_start["c_fol"],
        "wood_ratio": out.pools_end["c_woo"] / out.pools_start["c_woo"],
        "soil_ratio": out.pools_end["c_som"] / out.pools_start["c_som"],
        "nee_cycle": out.nee @ weights,   # (draws, 365)
        "foliage_rs": steady_state_ratio(out.pools_in["c_fol"], out.pools_out["c_fol"],
                                         out.pools_yearly["c_fol"].mean(axis=-1),
                                         out.pools_start["c_fol"]),
    }


def run_all(params, drivers, latitude: float, frost: float | None, years, weights,
            batch_size: int) -> dict[str, np.ndarray]:
    """Run every draw in batches and stack the per-draw summaries."""
    p = params.as_arrays()
    n = p["f_auto"].shape[0]
    parts = []
    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        batch = type(params)(**{k: v[start:stop] for k, v in p.items()})
        with np.errstate(all="ignore"):   # impossible draws may overflow; they are counted
            out = run(batch, drivers, latitude_deg=latitude, frost_threshold_degc=frost)
            parts.append(summarise_batch(out, years, weights))
        print(f"  draws {start + 1}-{stop} of {n} done")
    return {k: np.concatenate([part[k] for part in parts]) for k in parts[0]}


# ---------------------------------------------------------------------------
# Applying D40
# ---------------------------------------------------------------------------

def evaluate(m: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """True/False per draw for every D40 check."""
    checks = {
        "not a number": ~m["finite"],
        "a pool below zero": m["negative_pool"],
        "LAI above 12": m["lai_max"] > 12,
        "daily GPP above 20": m["gpp_daily_max"] > 20,
        "annual GPP above 3,000": m["gpp_annual_max"] > 3000,
    }
    tier1_fail = np.any(list(checks.values()), axis=0)
    tier2 = {name: (m[name] >= low) & (m[name] <= high) & ~tier1_fail
             for name, (low, high, _) in TIER2_BANDS.items()}
    return {"tier1_reasons": checks, "tier1_fail": tier1_fail, "tier2": tier2,
            "tier2_all": np.all(list(tier2.values()), axis=0)}


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


def prior_scale(name: str, values: np.ndarray) -> np.ndarray:
    """A sampled quantity on the scale its prior lives on (log, logit or plain)."""
    kind = PRIORS[name][0] if name in PRIORS else "dirichlet"
    if kind == "lognormal":
        return np.log(values)
    if kind in ("beta", "dirichlet"):
        v = np.clip(values, 1e-12, 1 - 1e-12)
        return np.log(v / (1 - v))
    return values


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
        x = prior_scale(name, values)
        links[name] = float((x[failed].mean() - x[~failed].mean()) / x.std())
    print(f"\n  All sampled quantities, linked to '{label}' (failing minus passing):")
    for name, diff in sorted(links.items(), key=lambda kv: -abs(kv[1])):
        direction = "higher" if diff > 0 else "lower"
        flag = "  <-- worth a look" if abs(diff) >= 0.2 else ""
        print(f"      {name:22s} {diff:+.2f}  (failing draws have {direction} values){flag}")
    return links


def parameter_table(raw, tier1_links, tier2_links, config) -> None:
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
    path = resolve_path(config["paths"]["figures"]) / "prior_predictive_parameters.csv"
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


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def figure_l(m: dict[str, np.ndarray], observed: pd.DataFrame, config: dict[str, Any]) -> None:
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
    ax.set_title(f"NEE through the year implied by the priors, {start}-{end}")
    ax.set_ylabel(f"NEE ({UNIT_FLUX})\nnegative = forest takes up carbon")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncols=2)
    save_figure(fig, "fig_L_prior_predictive_nee", config)


def figure_m(m: dict[str, np.ndarray], config: dict[str, Any]) -> None:
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
    fig.suptitle("What the priors imply about the forest (prior predictive, D40)", fontsize=10)
    fig.tight_layout()
    save_figure(fig, "fig_M_prior_predictive_checks", config)


# ---------------------------------------------------------------------------

def main() -> None:
    config = load_config()
    apply_style()
    settings = config["prior_predictive"]

    data = load_model_block(config, block="calibration")
    drivers = drivers_from_frame(data, config)
    years = data.index.year.to_numpy()
    weights = doy_weights(pd.DatetimeIndex(data.index))
    latitude = float(config["site"]["latitude_deg"])

    rng = np.random.default_rng(int(settings["seed"]))
    params, raw = sample_prior(int(settings["n_draws"]), rng)

    heading(f"Running {settings['n_draws']} prior draws over {drivers.n_days} days "
            "(main model: no frost cut-off)")
    main_run = run_all(params, drivers, latitude, None, years, weights, int(settings["batch_size"]))
    main_result = evaluate(main_run)
    report(main_result, "Prior predictive check, main model (D40)")
    tier1_links = drivers_of_failure(raw, main_result["tier1_fail"], "Tier 1 failure")
    tier2_links = drivers_of_failure(raw, ~main_result["tier2_all"], "Tier 2 failure")
    parameter_table(raw, tier1_links, tier2_links, config)

    heading("Same draws with the frost switch on (D24), reported alongside")
    frost = float(config["model"]["frost_cutoff"]["threshold_degc"])
    frost_run = run_all(params, drivers, latitude, frost, years, weights, int(settings["batch_size"]))
    report(evaluate(frost_run), f"Prior predictive check, frost cut-off at {frost} degC")

    weighted_report(main_result, edc_weights(main_run),
                    "With the soft foliage constraint (D44), main model, by re-weighting")
    weighted_report(evaluate(frost_run), edc_weights(frost_run),
                    "With the soft foliage constraint (D44), frost switch on, by re-weighting")

    heading("Medians of the plausible draws (main model)")
    ok = ~main_result["tier1_fail"]
    for name in TIER2_BANDS:
        print(f"  {name:14s} median {np.median(main_run[name][ok]):9.2f}")

    qc_pass = data[MASK_NEE] & data[MASK_QC]
    observed = climatology(data[config["columns"]["nee"]].where(qc_pass))
    figure_l(main_run, observed, config)
    figure_m(main_run, config)


if __name__ == "__main__":
    main()