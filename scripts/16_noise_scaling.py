"""RQ2 noise-scaling experiment (D71): what happens to the posterior when RANDUNC is
multiplied by 0.5 and by 2? Reads the tables and fits that scripts 04 and 05 saved; nothing is sampled.

For every quantity and each arm:
    width ratio at the three levels c = 0.5, 1, 2   (posterior 90% width / prior 90% width, prior scale, D32)
    response = slope of log(width ratio) against log(c)   (1: noise passes straight through; 0: it does not reach it)
Arm 1 (fake data) also has the prediction written down in D71:
    predicted response = 1 - r^2, with r the width ratio at c = 1;  it "holds" if |measured - predicted| <= 0.2.
Arm 2 (real data) also reports sigma_m at each level and the weight of a July day against a January day,
    (c^2 R_Jan^2 + sigma_m^2) / (c^2 R_Jul^2 + sigma_m^2),  with R_Jan = 0.06 and R_Jul = 0.27 (F43).

A run is used only if it passes the computation checks of D49 (no divergences, R-hat, ESS); a run that fails
is listed as failing and left out. Pool widths (end of 2010) are NOT in this script: they need a forward run.

Run from the repository root:
    python scripts/16_noise_scaling.py --fake fake_data_check_cold_MC_dense_ta995_d2000 \\
                                       --real fit_model_error_MC_cold_dense_ta995_d2000
(the _ru0p5 and _ru2 runs are found from these two level-1 names)

Saves in figures/:
    noise_scaling_fake.csv, noise_scaling_real.csv, fig_AF_noise_scaling_fake, fig_AG_noise_scaling_real
"""

import argparse
import importlib.util
from pathlib import Path

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.plotting import BLUE, GREY, GREY_LIGHT, INK_SOFT, ORANGE, apply_style, save_figure
from dalec.posterior import computation_checks, randunc_scale_tag, sampler_health

_here = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("talk_figures", _here / "talk_figures.py")
talk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(talk)
_spec15 = importlib.util.spec_from_file_location("trace_overview", _here / "15_trace_overview.py")
labels15 = importlib.util.module_from_spec(_spec15)
_spec15.loader.exec_module(labels15)

LEVELS = (0.5, 1.0, 2.0)
DERIVED = ("f_woo",)             # computed from the sampled fractions, not sampled
HOLDS_WITHIN = 0.2               # D71
R_JAN, R_JUL = 0.06, 0.27        # typical RANDUNC in January and July (F43), g C m-2 d-1


def level_label(level1: str, c: float) -> str:
    """File label of one level: the level-1 label with the _ru tag added at the end."""
    return level1 + randunc_scale_tag(c)


def table_path(label: str, figures: Path) -> Path:
    """Where script 04 (fake data) or 05 (real data) saved the table of a run."""
    if label.startswith("fake_data_check"):
        return figures / (label.replace("fake_data_check", "fake_data_recovery", 1) + ".csv")
    return figures / f"{label}_parameters.csv"


def response(widths: dict[float, float]) -> float:
    """Slope of log(width ratio) against log(c) over the levels."""
    c = np.log(list(widths))
    w = np.log(list(widths.values()))
    return float(np.polyfit(c, w, 1)[0])


def noise_scaling_table(tables: dict[float, pd.DataFrame], with_prediction: bool) -> pd.DataFrame:
    """One row per sampled quantity: width ratio at each level, measured response, and (arm 1) the prediction."""
    base = tables[1.0]
    names = [n for n in base.index if base.loc[n, "narrowing"] != "fixed" and n not in DERIVED]
    rows = []
    for n in names:
        widths = {c: float(tables[c].loc[n, "width_ratio"]) for c in tables}
        row = {"quantity": n, **{f"ratio_c{c:g}": w for c, w in widths.items()}, "response": response(widths)}
        if with_prediction:
            r = widths[1.0]
            row["predicted"] = 1.0 - r ** 2
            row["holds"] = abs(row["response"] - row["predicted"]) <= HOLDS_WITHIN
        rows.append(row)
    return pd.DataFrame(rows).set_index("quantity")


def day_weight(sigma_m: float, c: float) -> float:
    """Weight of a July day relative to a January day: (c^2 R_Jan^2 + s^2) / (c^2 R_Jul^2 + s^2)."""
    return (c ** 2 * R_JAN ** 2 + sigma_m ** 2) / (c ** 2 * R_JUL ** 2 + sigma_m ** 2)


def fake_figure(table: pd.DataFrame) -> plt.Figure:
    """Measured response against the prediction 1 - r^2, one point per quantity."""
    fig, ax = plt.subplots(figsize=(talk.HALF_WIDTH + 0.6, 5.6))
    lim = (-0.15, 1.1)
    xs = np.linspace(*lim, 10)
    ax.fill_between(xs, xs - HOLDS_WITHIN, xs + HOLDS_WITHIN, color=GREY_LIGHT, alpha=0.5, lw=0, zorder=0)
    ax.plot(xs, xs, color=GREY, lw=1, ls="--", zorder=1)
    colours = [BLUE if h else ORANGE for h in table["holds"]]
    ax.scatter(table["predicted"], table["response"], s=48, c=colours, edgecolors="white", linewidths=1, zorder=3)
    for n, row in table.iterrows():
        ax.annotate(labels15.thesis_label(n), (row["predicted"], row["response"]), xytext=(4, 3),
                    textcoords="offset points", fontsize=9, color=INK_SOFT)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(r"predicted response, $1 - r^2$")
    ax.set_ylabel("measured response (slope of log width on log noise)")
    ax.grid(False)
    fig.tight_layout()
    return fig


def real_figure(table: pd.DataFrame) -> plt.Figure:
    """Measured response per quantity, real data (sigma_m free)."""
    t = table["response"].sort_values()
    fig, ax = plt.subplots(figsize=(talk.HALF_WIDTH + 1.4, 0.34 * len(t) + 1.6))
    rows = range(len(t))
    ax.barh(rows, t, height=0.62, color=BLUE, zorder=2)
    ax.set_yticks(rows, [labels15.thesis_label(n) for n in t.index])
    ax.axvline(0, color=GREY, lw=1)
    ax.axvline(1, color=GREY, lw=1, ls="--")
    ax.set_xlabel("response to the noise (0: none, 1: passes straight through)")
    ax.grid(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    return fig


def run_arm(level1: str, config: dict, figures: Path, results: Path, with_prediction: bool):
    """Check D49 on all three levels, then build the table from those that pass.
    Returns the table (or None) and the posterior median of sigma_m at each level (real data only)."""
    tables, failed, sigma_m = {}, [], {}
    for c in LEVELS:
        label = level_label(level1, c)
        path = table_path(label, figures)
        nc = results / f"{label}.nc"
        if not path.exists() or not nc.exists():
            print(f"  c = {c:g}: MISSING ({path.name} or {nc.name})")
            failed.append(c)
            continue
        table = pd.read_csv(path, index_col=0)
        idata = az.from_netcdf(nc)
        health = sampler_health(idata)
        checks = computation_checks(table, health, config["diagnostics"])
        ok = all(checks.values())
        print(f"  c = {c:g}: {label}  D49 {'PASS' if ok else 'FAIL'}  (divergences {health['divergences']}, "
              f"max R-hat {table['r_hat'].max():.4f}, min ESS {min(table['ess_bulk'].min(), table['ess_tail'].min()):.0f})")
        if ok:
            tables[c] = table
            if "sigma_m" in idata.posterior.data_vars:      # not in the parameter table; read from the fit itself
                sigma_m[c] = float(np.median(np.asarray(idata.posterior["sigma_m"].values)))
        else:
            failed.append(c)
    if len(tables) < len(LEVELS):
        print(f"  not all three levels are usable (failed or missing: {failed}); nothing is interpreted for this arm")
        return None, sigma_m
    return noise_scaling_table(tables, with_prediction), sigma_m


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fake", help="level-1 fake-data check, without .nc")
    parser.add_argument("--real", help="level-1 real-data fit, without .nc")
    args = parser.parse_args()
    config = load_config()
    apply_style()
    plt.rcParams.update(talk.TALK_STYLE)
    figures = resolve_path(config["paths"]["figures"])
    results = resolve_path(config["paths"]["results"])

    if args.fake:
        print("Arm 1, fake data")
        table, _ = run_arm(args.fake, config, figures, results, with_prediction=True)
        if table is not None:
            print(table.round(3).to_string())
            print(f"  prediction holds (within {HOLDS_WITHIN}) for {int(table['holds'].sum())} of {len(table)} quantities")
            table.to_csv(figures / "noise_scaling_fake.csv", float_format="%.5g")
            save_figure(fake_figure(table), "fig_AF_noise_scaling_fake", config)
    if args.real:
        print("Arm 2, real data")
        table, sig = run_arm(args.real, config, figures, results, with_prediction=False)
        if table is not None:
            print(table.round(3).to_string())
            for c in LEVELS:
                print(f"  c = {c:g}: sigma_m {sig[c]:.3f}   July day weight against January {day_weight(sig[c], c):.2f}")
            table.to_csv(figures / "noise_scaling_real.csv", float_format="%.5g")
            save_figure(real_figure(table), "fig_AG_noise_scaling_real", config)


if __name__ == "__main__":
    main()
