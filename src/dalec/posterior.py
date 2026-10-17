"""Shared steps for every fit (fake data, script 04; real data, script 05).

Finding chain starts, summarising the posterior against the prior, checking the
sampler's health, and drawing the prior-versus-posterior figure. One definition,
so the fake-data check and the real fits are judged the same way.
"""

import importlib.util
from typing import Any

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.inference import SHARES
from dalec.model_numpy import Drivers, Parameters
from dalec.plausibility import evaluate, run_and_summarise
from dalec.plotting import BLUE, GREY_LIGHT, INK_SOFT, ORANGE, save_figure
from dalec.priors import (
    ACM_PRIORS, COLD_PRIORS, priors_for, priors_from_attrs, sample_acm, sample_cold, sample_prior,
    to_prior_scale,
)

# The 24 sampled quantities, in the order of the parameter groups (D34 to D39).
QUANTITIES = ("f_auto", *SHARES, "theta_roo", "theta_woo", "theta_lit", "theta_som", "h",
              "temperature_exponent", "d_onset", "cr_onset", "d_fall", "cr_fall", "c_lf",
              "lma", "ceff", "c_lab_0", "c_fol_0", "c_roo_0", "c_woo_0", "c_lit_0", "c_som_0")

# Words for how much the data narrowed each prior (D49; descriptive, not pass/fail).
NARROWING = ((0.5, "strongly informed"), (0.9, "partly informed"), (np.inf, "not informed"))


# ---------------------------------------------------------------------------
# Before sampling
# ---------------------------------------------------------------------------

def check_can_save() -> None:
    """Fail now, not after hours of sampling, if results cannot be written.
    Saving the posterior as NetCDF needs the netCDF4 package (environment.yml)."""
    if importlib.util.find_spec("netCDF4") is None:
        raise ImportError("saving the posterior needs the netcdf4 package: add it to "
                          "environment.yml and run conda env update")


def first_days(drivers: Drivers, n: int) -> Drivers:
    """The same drivers cut to the first n days (quick runs)."""
    return Drivers(t=drivers.t[:n], doy=drivers.doy[:n], t_max=drivers.t_max[:n],
                   t_min=drivers.t_min[:n], sw_in=drivers.sw_in[:n], co2=drivers.co2[:n])


def search_prior(drivers, years, latitude, frost, *, seed: int, n_candidates: int,
                 batch_size: int, cold: bool = False, free_a7: bool = False):
    """Draw prior parameter sets, run each over the block, and apply the D40 checks.

    Returns (params, raw, summary, checks). Callers pick from these: physically
    possible draws (Tier 1 passed) as chain starts (D45), and for the fake-data
    check, a fully plausible draw as the truth (D48). With `cold`, the draws come
    from the cold model's priors and include the two thresholds (D55, D57); with
    `free_a7`, they include the ACM light-use efficiency (D58).
    """
    rng = np.random.default_rng(seed)
    params, raw = sample_prior(n_candidates, rng, priors_for(cold))
    thresholds = None
    if cold:
        raw.update(sample_cold(n_candidates, rng))
        thresholds = tuple(raw[name] for name in COLD_PRIORS)
    a7 = None
    if free_a7:
        raw.update(sample_acm(n_candidates, rng))
        a7 = raw["acm_a7"]
    summary = run_and_summarise(params, drivers, latitude, frost, years, None,
                                batch_size, verbose=False, cold=thresholds, a7=a7)
    return params, raw, summary, evaluate(summary)


def possible_draws(checks: dict[str, Any], n: int, exclude: tuple[int, ...] = ()) -> list[int]:
    """The first n prior draws that pass Tier 1 (chain starts, D45)."""
    starts = [int(i) for i in np.flatnonzero(~checks["tier1_fail"]) if i not in exclude][:n]
    if len(starts) < n:
        raise RuntimeError("not enough physically possible prior draws for the chain starts; "
                           "raise n_candidates in the config")
    return starts


# ---------------------------------------------------------------------------
# After sampling
# ---------------------------------------------------------------------------

def sampler_health(idata) -> dict[str, float]:
    """Divergences, tree-depth limits, average trajectory length and step size."""
    stats = idata.sample_stats
    diverging = np.asarray(stats["diverging"].values)
    return {
        "divergences": int(diverging.sum()),
        "n_draws": int(diverging.size),
        "max_depth_hits": int(np.asarray(stats["reached_max_treedepth"].values).sum()),
        "mean_steps": float(np.asarray(stats["n_steps"].values).mean()),
        "mean_step_size": float(np.asarray(stats["step_size"].values).mean()),
    }


def parameter_table(idata, interval: float, seed: int,
                    truth: dict[str, float] | None = None) -> pd.DataFrame:
    """Prior and posterior interval for every quantity, how much the data
    narrowed the prior (widths compared on the prior's own scale, D32), and
    the convergence numbers. With `truth` (fake data), also whether the true
    value lies inside the posterior interval."""
    lo, hi = 50 * (1 - interval), 50 * (1 + interval)
    cold = "cold_t_shutdown" in idata.posterior.data_vars
    rng = np.random.default_rng(seed)
    _, prior = sample_prior(20_000, rng, priors_from_attrs(idata.attrs))   # the fit's own priors
    quantities = QUANTITIES
    if cold:
        prior.update(sample_cold(20_000, rng))
        quantities = quantities + tuple(COLD_PRIORS)
    if "acm_a7" in idata.posterior.data_vars:          # D58 diagnostic
        prior.update(sample_acm(20_000, rng))
        quantities = quantities + tuple(ACM_PRIORS)
    diagnostics = az.summary(idata, var_names=[q for q in quantities if
                                               float(np.ptp(np.asarray(idata.posterior[q].values))) > 0])
    rows = []
    for name in quantities:
        post = np.asarray(idata.posterior[name].values, dtype=float).reshape(-1)
        p_lo, p_mid, p_hi = np.percentile(prior[name], [lo, 50, hi])
        q_lo, q_mid, q_hi = np.percentile(post, [lo, 50, hi])
        scaled = {k: float(to_prior_scale(name, v)) for k, v in
                  {"p_lo": p_lo, "p_hi": p_hi, "q_lo": q_lo, "q_hi": q_hi}.items()}
        ratio = (scaled["q_hi"] - scaled["q_lo"]) / (scaled["p_hi"] - scaled["p_lo"])
        row = {"quantity": name}
        if truth is not None:
            row["truth"] = float(truth[name])
        row.update({
            "prior_low": p_lo, "prior_median": p_mid, "prior_high": p_hi,
            "posterior_low": q_lo, "posterior_median": q_mid, "posterior_high": q_hi,
        })
        if truth is not None:
            row["truth_inside"] = bool(q_lo <= truth[name] <= q_hi)
        row.update({
            "width_ratio": ratio,
            "narrowing": ("fixed" if q_hi == q_lo else
                          next(label for limit, label in NARROWING if ratio < limit)),
            "r_hat": float(diagnostics.loc[name, "r_hat"]) if name in diagnostics.index else np.nan,
            "ess_bulk": float(diagnostics.loc[name, "ess_bulk"]) if name in diagnostics.index else np.nan,
            "ess_tail": float(diagnostics.loc[name, "ess_tail"]) if name in diagnostics.index else np.nan,
        })
        rows.append(row)
    return pd.DataFrame(rows).set_index("quantity")


def print_table_and_health(table: pd.DataFrame, health: dict[str, float], heading,
                           title: str) -> None:
    heading(title)
    columns = [c for c in ("truth", "posterior_low", "posterior_median", "posterior_high",
                           "truth_inside", "width_ratio", "narrowing", "r_hat", "ess_bulk",
                           "ess_tail") if c in table.columns]
    with pd.option_context("display.width", 160, "display.max_columns", 20):
        print(table[columns].to_string(float_format=lambda v: f"{v:.4g}"))

    heading("Sampler health")
    print(f"  divergences: {health['divergences']} of {health['n_draws']} kept draws")
    print(f"  draws that hit the maximum tree depth: {health['max_depth_hits']}")
    print(f"  average gradient evaluations per draw: {health['mean_steps']:.0f} "
          f"(average step size {health['mean_step_size']:.4f})")
    print(f"  largest R-hat: {table['r_hat'].max():.4f} ({table['r_hat'].idxmax()})")
    print(f"  smallest effective sample size: bulk {table['ess_bulk'].min():.0f}, "
          f"tail {table['ess_tail'].min():.0f}")

    heading("How much NEE narrowed each prior (posterior width / prior width)")
    for label in [label for _, label in NARROWING] + ["fixed"]:
        names = table.index[table["narrowing"] == label].tolist()
        if names or label != "fixed":
            print(f"  {label:18s} {len(names):2d}: {', '.join(names)}")


def computation_checks(table: pd.DataFrame, health: dict[str, float],
                       diagnostics: dict[str, Any]) -> dict[str, bool]:
    """D49 criteria 1 to 3: the sampler must pass these before any result is read."""
    return {
        "no divergences": health["divergences"] == 0,
        f"R-hat below {diagnostics['max_rhat']} for every quantity":
            table["r_hat"].max() < diagnostics["max_rhat"],
        f"effective sample size at least {diagnostics['min_ess']} (bulk and tail)":
            min(table["ess_bulk"].min(), table["ess_tail"].min()) >= diagnostics["min_ess"],
    }


def interval_figure(table: pd.DataFrame, config: dict[str, Any], name: str, title: str) -> None:
    """Each quantity on its prior scale, shifted and stretched so that its prior
    interval runs from -1 to +1: grey = prior, blue = posterior, orange = truth
    (fake data only)."""
    has_truth = "truth" in table.columns
    fig, ax = plt.subplots(figsize=(7, 8))
    for row, (quantity, r) in enumerate(table.iloc[::-1].iterrows()):
        p_lo, p_hi = (float(to_prior_scale(quantity, r[k])) for k in ("prior_low", "prior_high"))
        centre, half = (p_lo + p_hi) / 2, (p_hi - p_lo) / 2
        def place(value):
            return (float(to_prior_scale(quantity, value)) - centre) / half
        ax.plot([-1, 1], [row, row], color=GREY_LIGHT, lw=7, solid_capstyle="butt", zorder=1)
        ax.plot([place(r["posterior_low"]), place(r["posterior_high"])], [row, row],
                color=BLUE, lw=3, solid_capstyle="butt", zorder=2)
        if has_truth:
            ax.plot(place(r["truth"]), row, "o", color=ORANGE, ms=5, zorder=3)
        else:
            ax.plot(place(r["posterior_median"]), row, "|", color=BLUE, ms=9, mew=2, zorder=3)
    ax.set_yticks(range(len(table)))
    ax.set_yticklabels(table.index[::-1], fontsize=8)
    ax.axvline(0, color=INK_SOFT, lw=0.6)
    ax.set_xlabel("position within the prior interval (-1 = low end, +1 = high end)")
    ax.set_title(title, fontsize=10)
    handles = [plt.Line2D([], [], color=GREY_LIGHT, lw=7, label="prior"),
               plt.Line2D([], [], color=BLUE, lw=3, label="posterior")]
    if has_truth:
        handles.append(plt.Line2D([], [], color=ORANGE, marker="o", lw=0, label="true value"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.08),
              ncols=len(handles))
    ax.grid(axis="y", visible=False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    save_figure(fig, name, config)


def posterior_parameters(idata, n: int, seed: int):
    """n posterior draws, chosen at random across all chains, for running the
    numpy model forward (posterior predictive checks).

    Returns (Parameters, sigma_m, cold, a7): sigma_m is the model error per draw
    (zeros without it); cold is (t_shutdown, t_full) per draw in the cold-limit
    model (D55), otherwise None; a7 is the sampled ACM light-use efficiency per
    draw (D58 diagnostic), otherwise None.
    """
    post = idata.posterior
    total = int(np.asarray(post["f_auto"].values).size)
    pick = np.sort(np.random.default_rng(seed).choice(total, size=min(n, total), replace=False))

    def draws(name):
        return np.asarray(post[name].values, dtype=float).reshape(-1)[pick]

    values = {f: draws(f) for f in Parameters.__dataclass_fields__}
    sigma_m = draws("sigma_m") if "sigma_m" in post.data_vars else np.zeros(pick.size)
    cold = ((draws("cold_t_shutdown"), draws("cold_t_full"))
            if "cold_t_shutdown" in post.data_vars else None)
    a7 = draws("acm_a7") if "acm_a7" in post.data_vars else None
    return Parameters(**values), sigma_m, cold, a7