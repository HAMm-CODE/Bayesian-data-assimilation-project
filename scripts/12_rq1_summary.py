"""RQ1 in one table and one figure: what NEE taught the chosen model about each quantity.

For every sampled quantity of a real-data fit, side by side:
    width_ratio_real   posterior 90% width / prior 90% width on the real data (1 = NEE taught nothing)
    width_ratio_fake   the same in the model's fake-data check(s), mean over the truths given
    position           where the posterior median sits in the prior, in units of the prior's
                       half 90% width, on the prior's own scale (log or logit, D32):
                       0 = at the prior median, +-1 = at the prior's 5% or 95% point
    outside_prior_90   |position| > 1
    own_prior, likelihood, reading   the own-prior sensitivity (D69 addition, F41)
Nothing is sampled; the script only reads the tables scripts 04, 05 and 11 saved.

The figure plots the fake-data width ratio against the real one, one point per
quantity, coloured and shaped by the own-prior reading. On the diagonal the real
data taught as much as the fake data did; above it, less; below it, more.

Run from the repository root, for example:
    python scripts/12_rq1_summary.py fit_model_error_MC_cold_dense_ta995_d2000 \
        --fake fake_data_check_cold_MC_dense_ta995_d2000 fake_data_check_cold_MC_truth2_dense_ta995_d2000

Saves (in figures/):
    rq1_<fit>.csv
    fig_AA_<fit>_fake_vs_real (pdf and png)
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.plotting import GREY, apply_style, save_figure
from dalec.priors import to_prior_scale

_spec = importlib.util.spec_from_file_location("talk_figures", Path(__file__).resolve().parent / "talk_figures.py")
talk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(talk)

# Three readings, three colours: the first three categorical slots of the dataviz
# reference palette, the only three that pass every pair for a scatter; each also
# has its own marker so colour is never the only cue.
READINGS = {
    "data decide; own prior hardly matters": ("data decide", "#2a78d6", "o"),
    "own prior and data pull against each other": ("prior and data pull against each other", "#eb6834", "s"),
    "own prior decides, data weak": ("prior decides, data weak", "#1baf7a", "^"),
    "neither moves it": ("neither", GREY, "x"),
}


# Label nudges (points) where two labels would overlap
LABEL_OFFSETS = {"cr_fall": (-36, -3), "c_lab_0": (6, -14), "d_fall": (6, 4), "f_auto": (6, -12)}


def position_in_prior(row: pd.Series) -> float:
    """Posterior median in units of the prior's half 90% width, on the prior's scale."""
    name = row.name
    lo, med, hi, post = (float(to_prior_scale(name, row[k]))
                         for k in ("prior_low", "prior_median", "prior_high", "posterior_median"))
    return (post - med) / ((hi - lo) / 2)


def rq1_table(fit_table: pd.DataFrame, fake_tables: list[pd.DataFrame], own: pd.DataFrame | None) -> pd.DataFrame:
    free = fit_table[fit_table["narrowing"] != "fixed"].copy()
    out = pd.DataFrame(index=free.index)
    out["width_ratio_real"] = free["width_ratio"].astype(float)
    if fake_tables:
        ratios = pd.concat([t["width_ratio"].astype(float) for t in fake_tables], axis=1)
        out["width_ratio_fake"] = ratios.reindex(out.index).mean(axis=1)
    out["position"] = free.apply(position_in_prior, axis=1)
    out["outside_prior_90"] = out["position"].abs() > 1
    if own is not None:
        own = own.rename(index=lambda n: n.removesuffix("_fraction"))   # cold thresholds are sampled as fractions
        out = out.join(own[["own_prior", "likelihood", "reading"]])
    return out


def figure(table: pd.DataFrame, label: str, n_truths: int, config) -> None:
    fig, ax = plt.subplots(figsize=(talk.HALF_WIDTH + 1.4, 6.2))
    ax.plot([0, 1], [0, 1], color=GREY, lw=1, ls="--", zorder=1)
    ax.text(0.97, 0.92, "real data taught less", ha="right", va="top", color="#52514e", fontsize=10.5)
    ax.text(0.97, 0.06, "real data taught more", ha="right", va="bottom", color="#52514e", fontsize=10.5)
    handles = {}
    for name, r in table.iterrows():
        text, colour, marker = READINGS.get(r.get("reading"), ("no sensitivity result", GREY, "o"))
        handles[text] = ax.scatter(r["width_ratio_fake"], r["width_ratio_real"], s=70, color=colour, marker=marker,
                                   edgecolors="white", linewidths=1.5, zorder=3, label=text)
        ax.annotate(talk.symbol_of(name), (r["width_ratio_fake"], r["width_ratio_real"]),
                    xytext=LABEL_OFFSETS.get(name, (5, 4)), textcoords="offset points", fontsize=10, color="#0b0b0b")
    order = [v[0] for v in READINGS.values()] + ["no sensitivity result"]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    truths = "one true forest" if n_truths == 1 else f"mean of {n_truths} true forests"
    ax.set_xlabel(f"width ratio, fake data ({truths})")
    ax.set_ylabel("width ratio, real data")
    shown = [k for k in order if k in handles]
    ax.legend([handles[k] for k in shown], shown, loc="upper left", frameon=False, fontsize=10)
    fig.tight_layout()
    save_figure(fig, f"fig_AA_{label}_fake_vs_real", config)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fit", help="real-data fit, without .nc")
    parser.add_argument("--fake", nargs="+", default=[], help="the same model's fake-data check(s), without .nc")
    args = parser.parse_args()
    config = load_config()
    apply_style()
    plt.rcParams.update(talk.TALK_STYLE)
    figures = resolve_path(config["paths"]["figures"])

    fit_table = pd.read_csv(figures / f"{args.fit}_parameters.csv", index_col=0)
    fakes = [pd.read_csv(figures / f"fake_data_recovery{run.removeprefix('fake_data_check')}.csv", index_col=0)
             for run in args.fake]
    own_path = figures / f"{args.fit}_prior_sensitivity_own.csv"
    own = pd.read_csv(own_path, index_col=0) if own_path.exists() else None
    table = rq1_table(fit_table, fakes, own)

    out = figures / f"rq1_{args.fit}.csv"
    table.to_csv(out, float_format="%.4g")
    with pd.option_context("display.width", 160):
        print(table.round(3).to_string())
    print(f"\n  saved {out}")
    if fakes:
        figure(table, args.fit, len(fakes), config)


if __name__ == "__main__":
    main()
