"""Fake-data check (Gelman et al. 2020, Section 4.1; DECISIONS.md D47 to D49).

The question: if we knew the right answer, would the fit find it?

1. Pick a "true" forest: the first prior draw that passes every D40 check
   (Tier 1 and Tier 2), so it is believable but not the prior centre (D48).
2. Run DALEC2 with it on the real 1997-2010 weather, and add random noise of
   the real size (RANDUNC) on the real assimilable days only. This is fake NEE
   whose true parameters we know.
3. Fit the full model (D47) to the fake NEE with NUTS, each chain starting from
   a different prior draw that passes Tier 1 (D45, D49).
4. Judge the fit by the criteria fixed in D49 before running: no divergences,
   R-hat below 1.01, effective sample size of at least 400, and the truth
   inside the 90% posterior interval for at least 20 of the 24 quantities.
   Also report how much each prior was narrowed by the data (RQ1 preview).

Run from the repository root:
    python scripts/04_fake_data_check.py --quick   # 2 years, 2 short chains: does it run? (minutes)
    python scripts/04_fake_data_check.py           # the real check (hours)
    python scripts/04_fake_data_check.py --mass-matrix dense   # D52 speed test
    python scripts/04_fake_data_check.py --cold    # the cold-limit model (D55, D59): truth,
                                                   # fake NEE and fit all include the cold limit
    python scripts/04_fake_data_check.py --cold --fix-set MB   # a reduced model (D63): the fixed
                                                   # quantities are fixed in the truth and in the fit
    python scripts/04_fake_data_check.py --cold --prior-set lit --fix-set MA   # literature priors (D64)

The quick run judges nothing; it only shows that every step works and prints
a rough time estimate for the full run.

Saves:
    results/fake_data_check[_cold][_<prior set>][_<fix set>][_<mass matrix>][_taNN][_quick].nc   all posterior draws (not tracked by git)
    figures/fake_data_recovery[...].csv   one row per quantity: truth, prior, posterior
    (the sampler tags appear only when --mass-matrix or --target-accept is given)
    N  fig_N_fake_data_recovery[...]      prior range, posterior range and truth
"""

import argparse
import time
from typing import Any

import numpy as np
from scipy import stats

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.inference import MASS_MATRIX_INIT, build_model, initial_values, sample
from dalec.model_numpy import Parameters, drivers_from_frame, run
from dalec.plotting import apply_style, heading
from dalec.priors import COLD_PRIORS, PRIOR_SETS, fixed_values, priors_for, priors_to_text
from dalec.posterior import (
    QUANTITIES, check_can_save, computation_checks, first_days, interval_figure,
    parameter_table, possible_draws, print_table_and_health, sampler_health, search_prior,
)


# ---------------------------------------------------------------------------
# 1 and 2: the true forest, the fake data and the starting points
# ---------------------------------------------------------------------------

def choose_truth_and_starts(drivers, years, latitude, frost, n_chains, config, cold=False, fixed=None,
                            prior_set=None):
    """Search prior draws: the first fully plausible one is the truth (D48);
    the next ones that are physically possible are the chain starts (D45)."""
    fake = config["fake_data"]
    params, raw, summary, checks = search_prior(
        drivers, years, latitude, frost, seed=int(fake["seed"]),
        n_candidates=int(fake["n_candidates"]),
        batch_size=int(config["prior_predictive"]["batch_size"]), cold=cold, fixed=fixed,
        prior_set=prior_set)
    plausible = np.flatnonzero(checks["tier2_all"])
    possible = np.flatnonzero(~checks["tier1_fail"])
    if plausible.size == 0:
        raise RuntimeError("no fully plausible prior draw found; raise fake_data.n_candidates")
    truth = int(plausible[0])
    starts = possible_draws(checks, n_chains, exclude=(truth,))
    print(f"  searched {params.as_arrays()['f_auto'].size} prior draws: "
          f"{possible.size} physically possible, {plausible.size} fully plausible")
    print(f"  true forest = draw {truth}; chain starts = draws {starts}")
    print("  the true forest, 1997-2010 averages: "
          f"GPP {summary['gpp'][truth]:.0f}, Rh {summary['rh'][truth]:.0f}, "
          f"NEE {summary['nee'][truth]:.0f} g C m-2 yr-1, peak LAI {summary['peak_lai'][truth]:.1f}, "
          f"foliage end/start {summary['foliage_ratio'][truth]:.2f}")
    if cold:
        print("  true cold limit: no GPP below "
              f"{raw['cold_t_shutdown'][truth]:.1f} degC, full GPP above {raw['cold_t_full'][truth]:.1f} degC")
    return params, raw, truth, starts


# ---------------------------------------------------------------------------
# 4: judging the fit
# ---------------------------------------------------------------------------

def allowed_misses(k: int, p_miss: float, tail: float) -> int:
    """Smallest m with P(more than m of k truths outside the interval) <= tail, when each
    truth falls outside with probability p_miss (D63). 4 for k = 22 to 26, 3 for k = 14."""
    m = 0
    while stats.binom.sf(m, k, p_miss) > tail:
        m += 1
    return m


def report(table, health, config: dict[str, Any], judged: bool) -> None:
    print_table_and_health(table, health, heading,
                           "Recovery of the true values (prior scale: log, logit or plain, D32)")
    if not judged:
        heading("Quick run: nothing is judged (too short); see the time estimate above")
        return
    fake = config["fake_data"]
    heading("D49 criteria, fixed before running")
    checks = computation_checks(table, health, config["diagnostics"])
    # Only the sampled quantities count; quantities held fixed (D63) are left out.
    # Misses allowed (D63): the smallest m with P(more than m misses) <= miss_tail for k
    # quantities at the 90% interval. This gives 4 for 22 to 26 quantities, the same as
    # D49 (20 of 24) and D59 (22 of 26), and 3 for the 14 of M0.
    free = table[table["narrowing"] != "fixed"]
    k = len(free)
    allowed = allowed_misses(k, 1.0 - float(fake["interval"]), float(fake["miss_tail"]))
    required = k - allowed
    inside = int(free["truth_inside"].sum())
    checks[f"truth inside the {fake['interval']:.0%} interval for at least "
           f"{required} of {k}"] = inside >= required
    for text, ok in checks.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {text}")
    print(f"  truth inside for {inside} of {k}; "
          f"outside: {', '.join(free.index[~free['truth_inside'].astype(bool)]) or 'none'}")
    print(f"\n  OVERALL: {'PASS' if all(checks.values()) else 'FAIL'}")


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true",
                        help="2 years and short chains: checks that everything runs")
    parser.add_argument("--mass-matrix", choices=tuple(MASS_MATRIX_INIT), default=None,
                        help="override sampler.mass_matrix from the config (D52 speed test)")
    parser.add_argument("--target-accept", type=float, default=None,
                        help="override sampler.target_accept (D53); higher = smaller, safer steps")
    parser.add_argument("--cold", action="store_true",
                        help="check the cold-limit model (D55, D59)")
    parser.add_argument("--fix-set", default=None,
                        help="a reduced model from config fix_sets (D63), e.g. --fix-set MB")
    parser.add_argument("--prior-set", choices=tuple(PRIOR_SETS), default=None,
                        help="use a named prior set (D64), e.g. --prior-set lit for the literature priors")
    args = parser.parse_args()

    check_can_save()
    config = load_config()
    apply_style()
    fake, settings = config["fake_data"], dict(config["sampler"])
    if args.quick:
        settings.update({k: v for k, v in settings["quick"].items() if k != "years"})
    mass_matrix = args.mass_matrix or settings["mass_matrix"]
    if args.target_accept is not None:
        settings["target_accept"] = args.target_accept
    # Output names carry a sampler tag only when the setting was given on the
    # command line; the config defaults (D53) carry none.
    ta_tag = "" if args.target_accept is None else f"_ta{round(100 * float(settings['target_accept']))}"
    mm_tag = "" if args.mass_matrix is None else f"_{args.mass_matrix}"
    priors = priors_for(args.cold, args.prior_set)
    prior_tag = f"_{args.prior_set}" if args.prior_set not in (None, "site") else ""
    fixed = None
    if args.fix_set:                                       # D63: a named reduced model
        if args.fix_set not in config["fix_sets"]:
            raise SystemExit(f"--fix-set {args.fix_set}: not in config fix_sets {sorted(config['fix_sets'])}")
        fixed = fixed_values(config["fix_sets"][args.fix_set], priors)
    suffix = (("_cold" if args.cold else "") + prior_tag + (f"_{args.fix_set}" if args.fix_set else "")
              + mm_tag + ta_tag + ("_quick" if args.quick else ""))

    data = load_model_block(config, block="calibration")
    drivers = drivers_from_frame(data, config)
    years = data.index.year.to_numpy()
    latitude = float(config["site"]["latitude_deg"])
    frost_cfg = config["model"]["frost_cutoff"]
    frost = float(frost_cfg["threshold_degc"]) if frost_cfg["enabled"] else None
    sigma = data[config["columns"]["sigma"]].to_numpy(dtype=float)
    use = data[MASK_USE].to_numpy(dtype=bool)

    heading("1. The true forest and the chain starts (D48, D45)")
    params, raw, truth, starts = choose_truth_and_starts(
        drivers, years, latitude, frost, int(settings["chains"]), config, cold=args.cold, fixed=fixed,
        prior_set=args.prior_set)
    if prior_tag:
        print(f"  prior set {args.prior_set!r} (D64), in the truth and in the fit")
    if fixed:
        print("  held fixed in the truth and in the fit (D63): " + ", ".join(sorted(fixed)))

    heading("2. Fake NEE: the true forest plus noise of the real size, on the real assimilable days")
    true_params = Parameters(**{k: v[truth] for k, v in params.as_arrays().items()})
    true_cold = tuple(float(raw[name][truth]) for name in COLD_PRIORS) if args.cold else None
    true_nee = run(true_params, drivers, latitude_deg=latitude, frost_threshold_degc=frost,
                   cold=true_cold).nee
    noise = np.random.default_rng(int(fake["seed"]) + 1).standard_normal(drivers.n_days)
    fake_nee = np.where(use, true_nee + np.where(use, sigma, 0.0) * noise, np.nan)
    n_days = drivers.n_days
    if args.quick:
        n_days = int(np.sum(years < years[0] + int(settings["quick"]["years"])))
        drivers = first_days(drivers, n_days)
    print(f"  {n_days} days, of which {int(use[:n_days].sum())} assimilable "
          f"(likelihood uses only these, D12)")

    heading("3. Building the model and sampling (D47, D49)")
    if config["likelihood"]["observation_model"] != "randunc_only":
        print("  note: the fake data have RANDUNC noise only, so the main observation model is used")
    model = build_model(drivers, fake_nee[:n_days], sigma[:n_days], use[:n_days],
                        latitude_deg=latitude, frost_threshold_degc=frost,
                        observation_model="randunc_only", cold_limit=args.cold, fixed=fixed,
                        prior_set=args.prior_set)
    print(f"  {settings['chains']} chains x ({settings['tune']} warm-up + {settings['draws']} kept), "
          f"target_accept {settings['target_accept']}, mass matrix {mass_matrix}")
    began = time.perf_counter()
    idata = sample(model, [initial_values(raw, i) for i in starts],
                   draws=int(settings["draws"]), tune=int(settings["tune"]),
                   target_accept=float(settings["target_accept"]), seed=int(settings["seed"]),
                   cores=int(settings["cores"]), mass_matrix=mass_matrix)
    elapsed = time.perf_counter() - began
    per_iteration = elapsed / (int(settings["tune"]) + int(settings["draws"]))
    print(f"  sampling took {elapsed / 60:.1f} min ({per_iteration:.2f} s per iteration per chain, "
          "chains side by side)")
    if args.quick:
        full = config["sampler"]
        scale = len(data) / n_days
        hours = per_iteration * scale * (int(full["tune"]) + int(full["draws"])) / 3600
        print(f"  rough estimate for the full run ({len(data)} days, {full['tune']} + {full['draws']} "
              f"per chain): about {hours:.1f} h, if the chains run side by side. Early warm-up is "
              "the slowest part, so treat this as a rough guide only.")

    out_dir = resolve_path(config["paths"]["results"])
    out_dir.mkdir(parents=True, exist_ok=True)
    idata.attrs["priors"] = priors_to_text(priors)   # D57, D64: which priors were used
    idata.to_netcdf(out_dir / f"fake_data_check{suffix}.nc")
    print(f"  saved {out_dir / f'fake_data_check{suffix}.nc'}")

    heading("4. Judging the fit (D49)")
    names = QUANTITIES + (tuple(COLD_PRIORS) if args.cold else ())
    truth_values = {name: float(raw[name][truth]) for name in names}
    table = parameter_table(idata, float(fake["interval"]), int(fake["seed"]) + 2, truth=truth_values)
    table_path = resolve_path(config["paths"]["figures"]) / f"fake_data_recovery{suffix}.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(table_path, float_format="%.6g")
    print(f"  saved {table_path}")
    report(table, sampler_health(idata), config, judged=not args.quick)
    interval_figure(table, config, f"fig_N_fake_data_recovery{suffix}",
                    f"Fake-data check{' (cold-limit model)' if args.cold else ''}: does the fit find the true values? "
                    f"({float(fake['interval']):.0%} intervals)")


if __name__ == "__main__":
    main()
