"""RQ1 in two figures: how much NEE narrowed each quantity, and whether the prior or the data decide it.

Nothing is sampled; the script only reads tables that scripts 05 and 11 saved.

Figure AB  uncertainty reduction, one bar per sampled quantity:
               reduction = (1 - width ratio) x 100
           where the width ratio is the posterior 90% width divided by the prior
           90% width, both on the prior's own scale (log or logit, D32). The dashed
           lines are the limits fixed before any run (D49): a width ratio of 0.9
           (10% reduction, "not informed" below) and 0.5 (50%, "strongly informed" above).
           Reads figures/<fit>_parameters.csv (script 05).

Figure AC  power-scaling sensitivity (Kallioinen et al. 2024; D69), two marks per quantity:
               its sensitivity to its OWN prior, and its sensitivity to the likelihood.
           A mark right of the dashed line (0.05, the ArviZ default) means "sensitive".
           Reads figures/<fit>_prior_sensitivity_own.csv (script 11 with --by-quantity).

Run from the repository root, for example:
    python scripts/13_rq1_figures.py fit_model_error_MC_cold_dense_ta995_d2000

Saves (in figures/, pdf and png):
    fig_AB_<fit>_uncertainty_reduction
    fig_AC_<fit>_prior_sensitivity      (only if the own-prior table exists)
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.plotting import BLUE, GREY, GREY_LIGHT, INK_SOFT, ORANGE, apply_style, save_figure

_spec = importlib.util.spec_from_file_location("talk_figures", Path(__file__).resolve().parent / "talk_figures.py")
talk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(talk)

LIMITS = ((10, "not informed below"), (50, "strongly informed above"))   # D49: width ratios 0.9 and 0.5
THRESHOLD = 0.05                                                          # power-scaling: "sensitive" above (ArviZ default)
DERIVED = ("f_woo",)                                                      # computed from the sampled fractions, not sampled


def reduction_table(table: pd.DataFrame) -> pd.Series:
    """Uncertainty reduction in percent for every sampled quantity, smallest first."""
    free = table[(table["narrowing"] != "fixed") & ~table.index.isin(DERIVED)]
    return (100 * (1 - free["width_ratio"].astype(float))).sort_values()


def reduction_figure(table: pd.DataFrame) -> plt.Figure:
    reduction = reduction_table(table)
    fig, ax = plt.subplots(figsize=(talk.HALF_WIDTH + 1.4, 0.34 * len(reduction) + 1.6))
    rows = range(len(reduction))
    ax.barh(rows, reduction, height=0.62, color=BLUE, zorder=2)
    ax.set_yticks(rows, [talk.symbol_of(n) for n in reduction.index])
    ax.set_ylim(-0.7, len(reduction) + 0.3)
    for x, text in LIMITS:
        ax.axvline(x, color=GREY, lw=1, ls="--", zorder=1)
        ax.text(x + 1, len(reduction) - 0.3, f"{x}%: {text}", color=INK_SOFT, fontsize=9.5, va="bottom")
    ax.axvline(0, color=GREY, lw=1, zorder=1)
    for y, r in zip(rows, reduction):
        ax.text(r + 1.2 if r >= 0 else r - 1.2, y, f"{round(r):d}%", va="center",
                ha="left" if r >= 0 else "right", color=INK_SOFT, fontsize=9.5)
    ax.set_xlim(min(-10.0, float(reduction.min()) - 8), 106)
    ax.set_xlabel("uncertainty reduction (%)")
    ax.grid(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    return fig


def sensitivity_table(own: pd.DataFrame) -> pd.DataFrame:
    """The two sensitivities per quantity, sorted by the sensitivity to the own prior."""
    own = own.rename(index=lambda n: str(n).removesuffix("_fraction"))   # cold thresholds are sampled as fractions
    return own[["own_prior", "likelihood"]].astype(float).sort_values("own_prior")


def sensitivity_figure(own: pd.DataFrame) -> plt.Figure:
    table = sensitivity_table(own)
    fig, ax = plt.subplots(figsize=(talk.HALF_WIDTH + 1.4, 0.34 * len(table) + 1.6))
    rows = range(len(table))
    for y, (p, l) in zip(rows, table.to_numpy()):
        ax.plot([p, l], [y, y], color=GREY_LIGHT, lw=2, zorder=1)
    ax.scatter(table["own_prior"], rows, s=62, color=ORANGE, marker="s", edgecolors="white", linewidths=1.2,
               zorder=3, label="sensitivity to its own prior")
    ax.scatter(table["likelihood"], rows, s=62, color=BLUE, marker="o", edgecolors="white", linewidths=1.2,
               zorder=3, label="sensitivity to the data (likelihood)")
    ax.axvline(THRESHOLD, color=GREY, lw=1, ls="--", zorder=0)
    ax.text(THRESHOLD * 1.04, len(table) - 0.3, f"{THRESHOLD:g}: sensitive above", color=INK_SOFT, fontsize=9.5,
            va="bottom")
    ax.set_yticks(rows, [talk.symbol_of(n) for n in table.index])
    ax.set_ylim(-0.7, len(table) + 0.3)
    ax.set_xscale("log")
    low = max(float(table.to_numpy().min()), 1e-4)
    ax.set_xlim(low / 1.6, float(table.to_numpy().max()) * 1.6)
    ax.set_xlabel("power-scaling sensitivity (log scale)")
    ax.grid(False)
    ax.tick_params(axis="y", length=0)
    fig.legend(loc="lower center", ncols=2, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fit", help="real-data fit, without .nc")
    args = parser.parse_args()
    config = load_config()
    apply_style()
    plt.rcParams.update(talk.TALK_STYLE)
    figures = resolve_path(config["paths"]["figures"])

    table = pd.read_csv(figures / f"{args.fit}_parameters.csv", index_col=0)
    print(reduction_table(table).round(1).to_string())
    save_figure(reduction_figure(table), f"fig_AB_{args.fit}_uncertainty_reduction", config)

    own_path = figures / f"{args.fit}_prior_sensitivity_own.csv"
    if own_path.exists():
        own = pd.read_csv(own_path, index_col=0)
        save_figure(sensitivity_figure(own), f"fig_AC_{args.fit}_prior_sensitivity", config)
    else:
        print(f"  no {own_path.name}: run script 11 with --by-quantity first; figure AC skipped")


if __name__ == "__main__":
    main()
