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

The quick run judges nothing; it only shows that every step works and prints
a rough time estimate for the full run.

Saves:
    results/fake_data_check[_quick].nc        all posterior draws (not tracked by git)
    figures/fake_data_recovery[_quick].csv    one row per quantity: truth, prior, posterior
    N  fig_N_fake_data_recovery[_quick]       prior range, posterior range and truth
"""

import argparse
import importlib.util
import time
from typing import Any

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.inference import SHARES, build_model, initial_values, sample
from dalec.model_numpy import Drivers, Parameters, drivers_from_frame, run
from dalec.plausibility import evaluate, run_and_summarise
from dalec.plotting import BLUE, GREY_LIGHT, INK_SOFT, ORANGE, apply_style, heading, save_figure
from dalec.priors import sample_prior, to_prior_scale

# The 24 sampled quantities, in the order of the parameter groups (D34 to D39).
QUANTITIES = ("f_auto", *SHARES, "theta_roo", "theta_woo", "theta_lit", "theta_som", "h",
              "temperature_exponent", "d_onset", "cr_onset", "d_fall", "cr_fall", "c_lf",
              "lma", "ceff", "c_lab_0", "c_fol_0", "c_roo_0", "c_woo_0", "c_lit_0", "c_som_0")

# Words for how much the data narrowed each prior (D49; descriptive, not pass/fail).
NARROWING = ((0.5, "strongly informed"), (0.9, "partly informed"), (np.inf, "not informed"))


# ---------------------------------------------------------------------------
# 1 and 2: the true forest, the fake data and the starting points
# ---------------------------------------------------------------------------

def choose_truth_and_starts(drivers, years, latitude, frost, n_chains, config):
    """Search prior draws: the first fully plausible one is the truth (D48);
    the next ones that are physically possible are the chain starts (D45)."""
    fake = config["fake_data"]
    rng = np.random.default_rng(int(fake["seed"]))
    params, raw = sample_prior(int(fake["n_candidates"]), rng)
    summary = run_and_summarise(params, drivers, latitude, frost, years, None,
                                int(config["prior_predictive"]["batch_size"]), verbose=False)
    checks = evaluate(summary)
    plausible = np.flatnonzero(checks["tier2_all"])
    possible = np.flatnonzero(~checks["tier1_fail"])
    if plausible.size == 0:
        raise RuntimeError("no fully plausible prior draw found; raise fake_data.n_candidates")
    truth = int(plausible[0])
    starts = [int(i) for i in possible if i != truth][:n_chains]
    if len(starts) < n_chains:
        raise RuntimeError("not enough physically possible prior draws for the chain starts; "
                           "raise fake_data.n_candidates")
    print(f"  searched {params.as_arrays()['f_auto'].size} prior draws: "
          f"{possible.size} physically possible, {plausible.size} fully plausible")
    print(f"  true forest = draw {truth}; chain starts = draws {starts}")
    print("  the true forest, 1997-2010 averages: "
          f"GPP {summary['gpp'][truth]:.0f}, Rh {summary['rh'][truth]:.0f}, "
          f"NEE {summary['nee'][truth]:.0f} g C m-2 yr-1, peak LAI {summary['peak_lai'][truth]:.1f}, "
          f"foliage end/start {summary['foliage_ratio'][truth]:.2f}")
    return params, raw, truth, starts


def check_can_save() -> None:
    """Fail now, not after hours of sampling, if results cannot be written.
    Saving the posterior as NetCDF needs the netCDF4 package (environment.yml)."""
    if importlib.util.find_spec("netCDF4") is None:
        raise ImportError("saving the posterior needs the netcdf4 package: add it to "
                          "environment.yml and run conda env update")


def first_days(drivers: Drivers, n: int) -> Drivers:
    """The same drivers cut to the first n days (quick run)."""
    return Drivers(t=drivers.t[:n], doy=drivers.doy[:n], t_max=drivers.t_max[:n],
                   t_min=drivers.t_min[:n], sw_in=drivers.sw_in[:n], co2=drivers.co2[:n])


# ---------------------------------------------------------------------------
# 4: judging the fit
# ---------------------------------------------------------------------------

def recovery_table(idata, raw, truth: int, interval: float, seed: int) -> pd.DataFrame:
    """Truth, prior and posterior for every quantity, on the prior's own scale (D32)."""
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    _, prior = sample_prior(20_000, np.random.default_rng(seed))
    diagnostics = az.summary(idata, var_names=list(QUANTITIES))
    rows = []
    for name in QUANTITIES:
        post = np.asarray(idata.posterior[name].values, dtype=float).reshape(-1)
        p_lo, p_mid, p_hi = np.percentile(prior[name], [lo, 50, hi])
        q_lo, q_mid, q_hi = np.percentile(post, [lo, 50, hi])
        true_value = float(raw[name][truth])
        scaled = {k: float(to_prior_scale(name, v)) for k, v in
                  {"p_lo": p_lo, "p_hi": p_hi, "q_lo": q_lo, "q_hi": q_hi}.items()}
        ratio = (scaled["q_hi"] - scaled["q_lo"]) / (scaled["p_hi"] - scaled["p_lo"])
        rows.append({
            "quantity": name, "truth": true_value,
            "prior_low": p_lo, "prior_median": p_mid, "prior_high": p_hi,
            "posterior_low": q_lo, "posterior_median": q_mid, "posterior_high": q_hi,
            "truth_inside": bool(q_lo <= true_value <= q_hi),
            "width_ratio": ratio,
            "narrowing": next(label for limit, label in NARROWING if ratio < limit),
            "r_hat": float(diagnostics.loc[name, "r_hat"]),
            "ess_bulk": float(diagnostics.loc[name, "ess_bulk"]),
            "ess_tail": float(diagnostics.loc[name, "ess_tail"]),
        })
    return pd.DataFrame(rows).set_index("quantity")


def report(table: pd.DataFrame, idata, config: dict[str, Any], judged: bool) -> None:
    fake = config["fake_data"]
    stats = idata.sample_stats
    divergences = int(np.asarray(stats["diverging"].values).sum())
    max_depth = int(np.asarray(stats["reached_max_treedepth"].values).sum())
    n_draws = int(np.asarray(stats["diverging"].values).size)

    heading("Recovery of the true values (prior scale: log, logit or plain, D32)")
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(table[["truth", "posterior_low", "posterior_median", "posterior_high",
                     "truth_inside", "width_ratio", "narrowing", "r_hat", "ess_bulk",
                     "ess_tail"]].to_string(float_format=lambda v: f"{v:.4g}"))

    heading("Sampler health")
    print(f"  divergences: {divergences} of {n_draws} kept draws")
    print(f"  draws that hit the maximum tree depth: {max_depth}")
    print(f"  largest R-hat: {table['r_hat'].max():.4f} ({table['r_hat'].idxmax()})")
    print(f"  smallest effective sample size: bulk {table['ess_bulk'].min():.0f}, "
          f"tail {table['ess_tail'].min():.0f}")

    heading("How much NEE narrowed each prior (posterior width / prior width)")
    for _, label in NARROWING:
        names = table.index[table["narrowing"] == label].tolist()
        print(f"  {label:18s} {len(names):2d}: {', '.join(names)}")

    if not judged:
        heading("Quick run: nothing is judged (too short); see the time estimate above")
        return
    heading("D49 criteria, fixed before running")
    checks = {
        "no divergences": divergences == 0,
        f"R-hat below {fake['max_rhat']} for every quantity": table["r_hat"].max() < fake["max_rhat"],
        f"effective sample size at least {fake['min_ess']} (bulk and tail)":
            min(table["ess_bulk"].min(), table["ess_tail"].min()) >= fake["min_ess"],
        f"truth inside the {fake['interval']:.0%} interval for at least "
        f"{fake['min_covered']} of {len(table)}": int(table["truth_inside"].sum()) >= fake["min_covered"],
    }
    for text, ok in checks.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {text}")
    print(f"  truth inside for {int(table['truth_inside'].sum())} of {len(table)}; "
          f"outside: {', '.join(table.index[~table['truth_inside']]) or 'none'}")
    print(f"\n  OVERALL: {'PASS' if all(checks.values()) else 'FAIL'}")


def figure_n(table: pd.DataFrame, config: dict[str, Any], suffix: str, interval: float) -> None:
    """Each quantity on its prior scale, shifted and stretched so that its prior
    interval runs from -1 to +1: grey = prior, blue = posterior, orange = truth."""
    fig, ax = plt.subplots(figsize=(7, 8))
    for row, (name, r) in enumerate(table.iloc[::-1].iterrows()):
        p_lo, p_hi = (float(to_prior_scale(name, r[k])) for k in ("prior_low", "prior_high"))
        centre, half = (p_lo + p_hi) / 2, (p_hi - p_lo) / 2
        def place(value):
            return (float(to_prior_scale(name, value)) - centre) / half
        ax.plot([-1, 1], [row, row], color=GREY_LIGHT, lw=7, solid_capstyle="butt", zorder=1)
        ax.plot([place(r["posterior_low"]), place(r["posterior_high"])], [row, row],
                color=BLUE, lw=3, solid_capstyle="butt", zorder=2)
        ax.plot(place(r["truth"]), row, "o", color=ORANGE, ms=5, zorder=3)
    ax.set_yticks(range(len(table)))
    ax.set_yticklabels(table.index[::-1], fontsize=8)
    ax.axvline(0, color=INK_SOFT, lw=0.6)
    ax.set_xlabel("position within the prior interval (-1 = low end, +1 = high end)")
    ax.set_title(f"Fake-data check: does the fit find the true values? ({interval:.0%} intervals)",
                 fontsize=10)
    handles = [plt.Line2D([], [], color=GREY_LIGHT, lw=7, label="prior"),
               plt.Line2D([], [], color=BLUE, lw=3, label="posterior"),
               plt.Line2D([], [], color=ORANGE, marker="o", lw=0, label="true value")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncols=3)
    ax.grid(axis="y", visible=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    save_figure(fig, f"fig_N_fake_data_recovery{suffix}", config)


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true",
                        help="2 years and short chains: checks that everything runs")
    args = parser.parse_args()

    check_can_save()
    config = load_config()
    apply_style()
    fake, settings = config["fake_data"], dict(config["sampler"])
    if args.quick:
        settings.update({k: fake["quick"][k] for k in ("chains", "tune", "draws")})
    suffix = "_quick" if args.quick else ""

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
        drivers, years, latitude, frost, int(settings["chains"]), config)

    heading("2. Fake NEE: the true forest plus noise of the real size, on the real assimilable days")
    true_params = Parameters(**{k: v[truth] for k, v in params.as_arrays().items()})
    true_nee = run(true_params, drivers, latitude_deg=latitude, frost_threshold_degc=frost).nee
    noise = np.random.default_rng(int(fake["seed"]) + 1).standard_normal(drivers.n_days)
    fake_nee = np.where(use, true_nee + np.where(use, sigma, 0.0) * noise, np.nan)
    n_days = drivers.n_days
    if args.quick:
        n_days = int(np.sum(years < years[0] + int(fake["quick"]["years"])))
        drivers = first_days(drivers, n_days)
    print(f"  {n_days} days, of which {int(use[:n_days].sum())} assimilable "
          f"(likelihood uses only these, D12)")

    heading("3. Building the model and sampling (D47, D49)")
    if config["likelihood"]["observation_model"] != "randunc_only":
        print("  note: the fake data have RANDUNC noise only, so the main observation model is used")
    model = build_model(drivers, fake_nee[:n_days], sigma[:n_days], use[:n_days],
                        latitude_deg=latitude, frost_threshold_degc=frost,
                        observation_model="randunc_only")
    print(f"  {settings['chains']} chains x ({settings['tune']} warm-up + {settings['draws']} kept), "
          f"target_accept {settings['target_accept']}")
    began = time.perf_counter()
    idata = sample(model, [initial_values(raw, i) for i in starts],
                   draws=int(settings["draws"]), tune=int(settings["tune"]),
                   target_accept=float(settings["target_accept"]), seed=int(settings["seed"]),
                   cores=int(settings["cores"]))
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
    idata.to_netcdf(out_dir / f"fake_data_check{suffix}.nc")
    print(f"  saved {out_dir / f'fake_data_check{suffix}.nc'}")

    heading("4. Judging the fit (D49)")
    table = recovery_table(idata, raw, truth, float(fake["interval"]), int(fake["seed"]) + 2)
    table_path = resolve_path(config["paths"]["figures"]) / f"fake_data_recovery{suffix}.csv"
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(table_path, float_format="%.6g")
    print(f"  saved {table_path}")
    report(table, idata, config, judged=not args.quick)
    figure_n(table, config, suffix, float(fake["interval"]))


if __name__ == "__main__":
    main()