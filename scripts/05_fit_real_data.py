"""Fit DALEC2 to the real FI-Hyy NEE, 1997-2010 (DECISIONS.md D51).

The same model (D47) and sampler settings (D49) that passed the fake-data check
(F16), now with the measured NEE_VUT_REF. Three parts:

1. Fit: NUTS, each chain starting from a different physically possible prior
   draw (D45). The posterior is saved.
2. Computation checks (D49 criteria 1 to 3: divergences, R-hat, effective
   sample size). If any fails, the results below must not be interpreted.
3. What the posterior says (Gelman et al. 2020, Sections 6.1 and 6.3), all
   descriptive, nothing tuned to it:
   - prior against posterior for every parameter (Figure O);
   - the forest the posterior implies: annual GPP, respiration, NEE, LAI and
     pool changes, next to the independent site values (RQ3);
   - posterior predictive checks: does the fitted model, with its error model,
     reproduce the observed NEE? Coverage of the 90% predictive interval,
     average misfit by season, the seasonal cycle (Figure P) and the monthly
     series (Figure Q).

Run from the repository root:
    python scripts/05_fit_real_data.py --quick        # 2 years, short chains: does it run? (minutes)
    python scripts/05_fit_real_data.py                # the main fit (about 30 h)
Sensitivity runs (planned runs table, Part 5) use the same script:
    --frost                                   frost switch on (D24)
    --cold                                    cold limit on photosynthesis, thresholds sampled (D55)
    --free-a7                                 ACM light-use efficiency estimated, a diagnostic (D58)
    --observation-model randunc_plus_model_error   model error added to sigma (L1)
    --mass-matrix dense                       sampler setting only (D52)
    --fix-set MB                              a reduced model from config fix_sets (D63)
    --prior-set lit                           the literature priors (D64)

Saves (the label is "fit" plus one tag per option, e.g. fit_frost):
    results/<label>.nc                 all posterior draws (not tracked by git)
    figures/<label>_parameters.csv     prior and posterior for every parameter
    figures/<label>_forest.csv         implied annual fluxes, LAI and pool changes
    O  fig_O_<label>_parameters        prior range against posterior range
    P  fig_P_<label>_seasonal_cycle    observed NEE through the year against the fitted model
    Q  fig_Q_<label>_monthly           observed monthly NEE against the fitted model, 1997-2010
"""

import argparse
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.inference import MASS_MATRIX_INIT, OBSERVATION_MODELS, build_model, initial_values, sample
from dalec.model_numpy import drivers_from_frame, run
from dalec.plausibility import TIER2_BANDS, evaluate, summarise_batch
from dalec.plotting import (
    BLUE, INK, INK_SOFT, UNIT_FLUX,
    apply_style, day_of_year_365, heading, save_figure, set_year_axis, smooth_cyclic,
)
from dalec.posterior import (
    check_can_save, computation_checks, interval_figure, parameter_table,
    possible_draws, posterior_parameters, print_table_and_health, sampler_health, search_prior,
)
from dalec.priors import PRIOR_SETS, fixed_values, prior_median, priors_for, priors_to_text

SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}

# Independent site values for context (Ilvesniemi et al. 2009; Kolari 2010; D40).
SITE_VALUES = {"gpp": "952 to 1,104", "rh": "290 to 370", "nee": "a sink",
               "peak_lai": "about 2.5 to 4"}


# ---------------------------------------------------------------------------
# Running the posterior forward
# ---------------------------------------------------------------------------

def posterior_runs(idata, drivers, latitude, frost, years, n_draws: int, seed: int,
                   batch_size: int):
    """Run the numpy model for n posterior draws. Returns daily NEE (draws x days),
    the extra model-error sd per draw (0 in the main model), and the D40 summaries."""
    params, sigma_m, cold, a7 = posterior_parameters(idata, n_draws, seed)
    arrays = params.as_arrays()
    nee, summaries = [], []
    for start in range(0, sigma_m.size, batch_size):
        stop = min(start + batch_size, sigma_m.size)
        batch = type(params)(**{k: v[start:stop] for k, v in arrays.items()})
        batch_cold = None if cold is None else (cold[0][start:stop], cold[1][start:stop])
        out = run(batch, drivers, latitude_deg=latitude, frost_threshold_degc=frost,
                  cold=batch_cold, a7=None if a7 is None else a7[start:stop])
        nee.append(out.nee)
        summary = summarise_batch(out, years)
        summary["ra"] = np.stack([out.ra[:, years == y].sum(axis=1)
                                  for y in np.unique(years)], axis=1).mean(axis=1)
        summaries.append(summary)
    summary = {k: np.concatenate([s[k] for s in summaries]) for k in summaries[0]}
    return np.concatenate(nee), sigma_m, summary


def predictive_draws(model_nee: np.ndarray, sigma: np.ndarray, sigma_m: np.ndarray,
                     use: np.ndarray, seed: int) -> np.ndarray:
    """Replicated observations: model NEE plus a fresh error of the size the
    likelihood assumes, on assimilable days (other days NaN)."""
    sd = np.sqrt(np.where(use, sigma, 0.0) ** 2 + sigma_m[:, None] ** 2)
    noise = np.random.default_rng(seed).standard_normal(model_nee.shape)
    return np.where(use, model_nee + sd * noise, np.nan)


# ---------------------------------------------------------------------------
# What the posterior says
# ---------------------------------------------------------------------------

def forest_table(summary: dict[str, np.ndarray], interval: float) -> pd.DataFrame:
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    names = {"gpp": "annual GPP", "ra": "annual autotrophic respiration",
             "rh": "annual heterotrophic respiration", "nee": "annual NEE (all days)",
             "peak_lai": "peak LAI", "foliage_ratio": "foliage pool, end / start",
             "wood_ratio": "wood pool, end / start", "soil_ratio": "soil pool, end / start"}
    rows = []
    for key, label in names.items():
        values = summary[key]
        rows.append({"quantity": label, "low": np.percentile(values, lo),
                     "median": np.median(values), "high": np.percentile(values, hi),
                     "site value": SITE_VALUES.get(key, "")})
    reco = summary["ra"] + summary["rh"]
    rows.insert(3, {"quantity": "annual ecosystem respiration", "low": np.percentile(reco, lo),
                    "median": np.median(reco), "high": np.percentile(reco, hi), "site value": ""})
    return pd.DataFrame(rows).set_index("quantity")


def predictive_checks(observed: np.ndarray, replicated: np.ndarray, model_nee: np.ndarray,
                      use: np.ndarray, months: np.ndarray, years: np.ndarray,
                      interval: float) -> None:
    """Coverage and average misfit, overall and by season; annual sums on assimilable days."""
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    inside = np.zeros(observed.size, dtype=bool)
    low, high = np.percentile(replicated[:, use], [lo, hi], axis=0)
    inside[use] = (observed[use] >= low) & (observed[use] <= high)
    residual = observed - model_nee.mean(axis=0)

    print(f"  share of assimilable days whose observed NEE lies inside the "
          f"{interval:.0%} predictive interval (expected about {interval:.0%} if the model "
          "and the error model are right):")
    print(f"      all seasons  {inside[use].mean():6.1%}   ({int(use.sum())} days)")
    for season, members in SEASONS.items():
        days = use & np.isin(months, members)
        print(f"      {season:12s} {inside[days].mean():6.1%}   ({int(days.sum())} days); "
              f"average observed minus model {residual[days].mean():+.2f} g C m-2 d-1")

    print("\n  annual NEE summed over assimilable days only (like for like), g C m-2 yr-1:")
    observed_years, model_years = [], []
    for y in np.unique(years):
        days = use & (years == y)
        observed_years.append(np.nansum(observed[days]))
        model_years.append(model_nee[:, days].sum(axis=1))
    model_years = np.stack(model_years, axis=1)        # draws x years
    for k, y in enumerate(np.unique(years)):
        print(f"      {y}: observed {observed_years[k]:7.0f}   model median "
              f"{np.median(model_years[:, k]):7.0f}  ({np.percentile(model_years[:, k], lo):.0f} "
              f"to {np.percentile(model_years[:, k], hi):.0f})")


def figure_p(observed, model_nee, use, index, config, name, interval) -> None:
    """Seasonal cycle on assimilable days: observed against the fitted model."""
    doy = day_of_year_365(pd.DatetimeIndex(index)).to_numpy()
    keep = use & ~np.isnan(doy)
    weights = np.zeros((len(index), 365))
    weights[np.flatnonzero(keep), doy[keep].astype(int) - 1] = 1.0
    counts = weights.sum(axis=0)
    weights = np.divide(weights, counts, out=np.zeros_like(weights), where=counts > 0)
    window = int(config["eda"]["smoothing_days"])

    def smooth(values):
        return smooth_cyclic(pd.Series(values, index=np.arange(1, 366)), window).to_numpy()

    obs_cycle = smooth(np.nan_to_num(observed) @ weights)
    model_cycles = np.stack([smooth(row) for row in model_nee @ weights])
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    days = np.arange(1, 366)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.fill_between(days, np.percentile(model_cycles, lo, axis=0),
                    np.percentile(model_cycles, hi, axis=0), color=BLUE, alpha=0.3, lw=0,
                    label=f"fitted model, {interval:.0%} interval of the posterior")
    ax.plot(days, np.median(model_cycles, axis=0), color=BLUE, lw=1.5, label="fitted model, median")
    ax.plot(days, obs_cycle, color=INK, lw=1.2, ls="--", label="observed")
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    set_year_axis(ax)
    first, last = pd.DatetimeIndex(index).year[[0, -1]]
    ax.set_title(f"NEE through the year, assimilable days {first}-{last} ({window}-day average)")
    ax.set_ylabel(f"NEE ({UNIT_FLUX})\nnegative = forest takes up carbon")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncols=3)
    save_figure(fig, name, config)


def figure_q(observed, model_nee, use, index, config, name, interval) -> None:
    """Monthly mean NEE on assimilable days: observed against the fitted model."""
    months = pd.DatetimeIndex(index).to_period("M")
    labels = months.unique()
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    obs_m, band_lo, band_mid, band_hi, when = [], [], [], [], []
    for m in labels:
        days = use & (months == m)
        if days.sum() < 5:          # too few measured days for a monthly mean
            continue
        model_mean = model_nee[:, days].mean(axis=1)
        obs_m.append(observed[days].mean())
        band_lo.append(np.percentile(model_mean, lo))
        band_mid.append(np.median(model_mean))
        band_hi.append(np.percentile(model_mean, hi))
        when.append(m.to_timestamp() + pd.Timedelta(days=14))
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.fill_between(when, band_lo, band_hi, color=BLUE, alpha=0.3, lw=0,
                    label=f"fitted model, {interval:.0%} interval")
    ax.plot(when, band_mid, color=BLUE, lw=1.2, label="fitted model, median")
    ax.plot(when, obs_m, "o", color=INK, ms=2.5, label="observed")
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    ax.set_title("Monthly mean NEE on assimilable days (months with at least 5 such days)")
    ax.set_ylabel(f"NEE ({UNIT_FLUX})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncols=3)
    save_figure(fig, name, config)


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true",
                        help="2 years and short chains: checks that everything runs")
    parser.add_argument("--frost", action="store_true", help="frost switch on (D24)")
    parser.add_argument("--cold", action="store_true",
                        help="add the cold limit on photosynthesis, two sampled thresholds (D55)")
    parser.add_argument("--free-a7", action="store_true",
                        help="estimate the ACM light-use efficiency a7 instead of fixing it (D58 diagnostic)")
    parser.add_argument("--observation-model", choices=OBSERVATION_MODELS, default=None,
                        help="override likelihood.observation_model from the config (L1)")
    parser.add_argument("--mass-matrix", choices=tuple(MASS_MATRIX_INIT), default=None,
                        help="override sampler.mass_matrix from the config (D52)")
    parser.add_argument("--fix", default="",
                        help="comma-separated parameters to hold at their prior median (D54), "
                             "e.g. --fix lma,c_lf")
    parser.add_argument("--target-accept", type=float, default=None,
                        help="override sampler.target_accept (D53); higher = smaller, safer steps")
    parser.add_argument("--fix-set", default=None,
                        help="a reduced model from config fix_sets (D63), e.g. --fix-set MB")
    parser.add_argument("--prior-set", choices=tuple(PRIOR_SETS), default=None,
                        help="use a named prior set (D64), e.g. --prior-set lit for the literature priors")
    args = parser.parse_args()

    check_can_save()
    config = load_config()
    apply_style()
    fit, settings, lik = config["fit"], dict(config["sampler"]), config["likelihood"]
    priors = priors_for(args.cold, args.prior_set)
    prior_tag = f"_{args.prior_set}" if args.prior_set not in (None, "site") else ""
    if args.quick:
        settings.update({k: v for k, v in settings["quick"].items() if k != "years"})
    observation_model = args.observation_model or lik["observation_model"]
    mass_matrix = args.mass_matrix or settings["mass_matrix"]
    frost_on = args.frost or bool(config["model"]["frost_cutoff"]["enabled"])
    frost = float(config["model"]["frost_cutoff"]["threshold_degc"]) if frost_on else None
    if args.target_accept is not None:
        settings["target_accept"] = args.target_accept
    # Output names carry a sampler tag only when the setting was given on the
    # command line; the config defaults (D53) carry none.
    ta_tag = "" if args.target_accept is None else f"_ta{round(100 * float(settings['target_accept']))}"
    mm_tag = "" if args.mass_matrix is None else f"_{args.mass_matrix}"
    fixed = {name: prior_median(name, priors) for name in filter(None, args.fix.split(","))}
    fix_tag = ("_fix-" + "-".join(sorted(fixed))) if fixed else ""
    if args.fix_set:                                       # D63: a named reduced model
        if args.fix_set not in config["fix_sets"]:
            raise SystemExit(f"--fix-set {args.fix_set}: not in config fix_sets {sorted(config['fix_sets'])}")
        fixed.update(fixed_values(config["fix_sets"][args.fix_set], priors))
        fix_tag += f"_{args.fix_set}"
    label = ("fit" + ("_model_error" if observation_model != "randunc_only" else "") + fix_tag
             + ("_frost" if frost_on else "") + ("_cold" if args.cold else "") + prior_tag
             + ("_a7" if args.free_a7 else "") + mm_tag + ta_tag
             + ("_quick" if args.quick else ""))
    interval = float(fit["interval"])
    batch_size = int(config["prior_predictive"]["batch_size"])

    data = load_model_block(config, block="calibration")
    if args.quick:
        years_all = data.index.year
        data = data[years_all < years_all[0] + int(settings["quick"]["years"])]
    drivers = drivers_from_frame(data, config)
    years = data.index.year.to_numpy()
    months = data.index.month.to_numpy()
    latitude = float(config["site"]["latitude_deg"])
    observed = data[config["columns"]["nee"]].to_numpy(dtype=float)
    sigma = data[config["columns"]["sigma"]].to_numpy(dtype=float)
    use = data[MASK_USE].to_numpy(dtype=bool)

    heading(f"1. Fitting the real NEE ({label})")
    print(f"  {drivers.n_days} days, of which {int(use.sum())} assimilable (D12)")
    print(f"  observation model {observation_model}; frost switch {'on' if frost_on else 'off'}; "
          f"mass matrix {mass_matrix}")
    if args.cold:
        print("  cold limit on photosynthesis on (D55): thresholds sampled")
    if args.free_a7:
        print("  ACM light-use efficiency a7 sampled, not fixed (D58 diagnostic)")
    if prior_tag:
        print(f"  prior set {args.prior_set!r} (D64): " + ", ".join(
            f"{k} {v}" for k, v in (PRIOR_SETS[args.prior_set] or {}).items()))
    if fixed:
        print("  held fixed at their prior medians (D54, D63): " + ", ".join(
            f"{k} = {np.round(v, 4).tolist() if np.ndim(v) else f'{v:.4g}'}" for k, v in fixed.items()))
    _, raw, _, checks = search_prior(drivers, years, latitude, frost, seed=int(fit["seed"]),
                                     n_candidates=int(fit["n_candidates"]), batch_size=batch_size,
                                     cold=args.cold, free_a7=args.free_a7, fixed=fixed,
                                     prior_set=args.prior_set)
    starts = possible_draws(checks, int(settings["chains"]))
    print(f"  chain starts = prior draws {starts} (each passes Tier 1, D45)")
    model = build_model(drivers, observed, sigma, use, latitude_deg=latitude,
                        frost_threshold_degc=frost, observation_model=observation_model,
                        model_error_prior_sd=float(lik["model_error_prior_sd"]), fixed=fixed,
                        cold_limit=args.cold, free_a7=args.free_a7, prior_set=args.prior_set)
    print(f"  {settings['chains']} chains x ({settings['tune']} warm-up + {settings['draws']} kept), "
          f"target_accept {settings['target_accept']}")
    began = time.perf_counter()
    idata = sample(model, [initial_values(raw, i) for i in starts],
                   draws=int(settings["draws"]), tune=int(settings["tune"]),
                   target_accept=float(settings["target_accept"]), seed=int(settings["seed"]),
                   cores=int(settings["cores"]), mass_matrix=mass_matrix)
    elapsed = time.perf_counter() - began
    print(f"  sampling took {elapsed / 3600:.1f} h "
          f"({elapsed / (int(settings['tune']) + int(settings['draws'])):.2f} s per iteration per chain)")
    idata.attrs["priors"] = priors_to_text(priors)   # D57, D64: which priors were used
    out_dir = resolve_path(config["paths"]["results"])
    out_dir.mkdir(parents=True, exist_ok=True)
    idata.to_netcdf(out_dir / f"{label}.nc")
    print(f"  saved {out_dir / f'{label}.nc'}")

    heading("2. Computation checks (D49 criteria 1 to 3)")
    table = parameter_table(idata, interval, int(fit["seed"]) + 2)
    health = sampler_health(idata)
    computation = computation_checks(table, health, config["diagnostics"])
    for text, ok in computation.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {text}")
    if args.quick:
        print("  (quick run: too short for these to mean anything)")
    elif not all(computation.values()):
        print("  -> at least one check FAILED: do not interpret the results below; "
              "find the cause first (D49, 'if it fails')")

    figures = resolve_path(config["paths"]["figures"])
    figures.mkdir(parents=True, exist_ok=True)
    table.to_csv(figures / f"{label}_parameters.csv", float_format="%.6g")
    print_table_and_health(table, health, heading, "3. Prior against posterior")
    interval_figure(table, config, f"fig_O_{label}_parameters",
                    f"What the real NEE says about each parameter ({interval:.0%} intervals)")

    heading("4. The forest the posterior implies (1997-2010 averages)")
    model_nee, sigma_m, summary = posterior_runs(
        idata, drivers, latitude, frost, years, int(fit["n_predictive_draws"]),
        int(fit["seed"]) + 3, batch_size)
    forest = forest_table(summary, interval)
    forest.to_csv(figures / f"{label}_forest.csv", float_format="%.6g")
    with pd.option_context("display.width", 160):
        print(forest.to_string(float_format=lambda v: f"{v:.3g}"))
    tier1 = evaluate(summary)["tier1_fail"].mean()
    print(f"\n  posterior draws that are physically impossible (D40 Tier 1): {tier1:.1%}")
    outside = {name: 1 - evaluate(summary)["tier2"][name].mean() for name in TIER2_BANDS}
    print("  share of posterior draws outside each D40 plausible band: "
          + ", ".join(f"{k} {v:.0%}" for k, v in outside.items()))

    heading("5. Posterior predictive checks (Gelman et al. 2020, Section 6.1)")
    replicated = predictive_draws(model_nee, sigma, sigma_m, use, int(fit["seed"]) + 4)
    predictive_checks(observed, replicated, model_nee, use, months, years, interval)
    figure_p(observed, model_nee, use, data.index, config, f"fig_P_{label}_seasonal_cycle", interval)
    figure_q(observed, model_nee, use, data.index, config, f"fig_Q_{label}_monthly", interval)


if __name__ == "__main__":
    main()
