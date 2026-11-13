"""RQ3 in two figures: can NEE find known values, and which quantities move together.

Nothing is sampled; the script reads saved fits and the recovery tables script 04 saved.

Figure AD  parameter recovery on fake data, one true forest: one small panel per
           sampled quantity with its prior (grey), its posterior (blue) and the
           true value (orange line). Panels are sorted by the width ratio, most
           narrowed first; a panel whose posterior 90% range misses the true
           value is marked "missed". Each panel shows its width ratio.
           Reads results/<run>.nc and figures/fake_data_recovery<...>.csv (script 04).
           With --intervals the older compact layout is drawn instead: one row
           per quantity on the prior's own scale, one panel per run given.

Figure AE  pairwise plot in one figure, to show correlations only: the quantities
           with at least one strong tie (|correlation| >= 0.5 on the prior scale,
           as in F44), at most --max of them. Lower triangle: one cloud of
           posterior draws per pair, with its correlation; a tinted panel is a
           strong tie. Diagonal: the posterior of each quantity. True values are
           never drawn here.
           Reads results/<fit>.nc.

Run from the repository root, for example:
    python scripts/14_rq3_figures.py \
        --recovery fake_data_check_cold_MC_dense_ta995_d2000 \
        --pairs fit_model_error_MC_cold_dense_ta995_d2000

Saves (in figures/, pdf and png):
    fig_AD_<run>_recovery
    fig_AE_<fit>_pairwise
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.plotting import BLUE, GREY, GREY_LIGHT, INK_SOFT, ORANGE, apply_style, save_figure
from dalec.priors import to_prior_scale

_spec = importlib.util.spec_from_file_location("talk_figures", Path(__file__).resolve().parent / "talk_figures.py")
talk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(talk)

STRONG = 0.5          # |correlation| at or above this is a strong tie (F44)
DERIVED = ("f_woo",)  # computed from the sampled fractions, not sampled
TINT = "#F1EAF8"      # light purple behind a strongly tied pair


# ---------------------------------------------------------------------------
# Figure AD: recovery
# ---------------------------------------------------------------------------

def recovery_density_figure(prior: dict[str, np.ndarray], posterior: dict[str, np.ndarray],
                            table: pd.DataFrame) -> plt.Figure:
    """Prior, posterior and true value of every sampled quantity of one fake-data check."""
    free = table[(table["narrowing"] != "fixed") & ~table.index.isin(DERIVED)]
    names = [n for n in free["width_ratio"].astype(float).sort_values().index if n in posterior]
    cols = 7 if len(names) > 20 else 5
    rows = int(np.ceil(len(names) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(talk.FULL_WIDTH, 1.3 * rows + 0.8), squeeze=False)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    for ax, name in zip(axes.flat, names):
        log_axis = talk.on_log_axis(name)
        if name in prior:
            x, y = talk.density(prior[name], log_axis)
            if x is not None:
                ax.fill_between(x, y, color=GREY_LIGHT, lw=0)
                ax.plot(x, y, color=GREY, lw=1)
        x, y = talk.density(posterior[name], log_axis)
        if x is not None:
            ax.plot(x, y, color=BLUE, lw=2)
        ax.axvline(float(free.loc[name, "truth"]), color=ORANGE, lw=2)
        if log_axis:
            ax.set_xscale("log")
        exponent = talk.scale_exponent(posterior[name]) if log_axis else None
        talk.few_ticks(ax, "x", 2, log_axis, exponent)
        ax.set_yticks([])
        ax.set_ylim(0, 1.32)
        ax.grid(False)
        ax.set_title(talk.with_scale(talk.short_label(name), exponent))
        text = f"{float(free.loc[name, 'width_ratio']):.2f}"
        if not bool(free.loc[name, "truth_inside"]):
            text += "  missed"
        ax.text(0.97, 0.96, text, transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color=INK_SOFT,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.2), zorder=5)
    handles = [plt.Line2D([], [], color=GREY, lw=7, alpha=0.4, label="prior"),
               plt.Line2D([], [], color=BLUE, lw=2.5, label="posterior"),
               plt.Line2D([], [], color=ORANGE, lw=2.5, label="true value"),
               plt.Line2D([], [], ls="", label="number: width ratio (posterior / prior 90% width)")]
    fig.legend(handles=handles, loc="lower center", ncols=4, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    return fig


def in_prior_units(table: pd.DataFrame) -> pd.DataFrame:
    """Truth and posterior of every sampled quantity in units of the prior's half 90%
    width, on the prior's scale: 0 = prior median, -1 and +1 = prior 5% and 95% points."""
    free = table[(table["narrowing"] != "fixed") & ~table.index.isin(DERIVED)]
    rows = {}
    for name, r in free.iterrows():
        lo, med, hi = (float(to_prior_scale(name, r[k])) for k in ("prior_low", "prior_median", "prior_high"))
        half = (hi - lo) / 2
        unit = lambda value: (float(to_prior_scale(name, value)) - med) / half
        rows[name] = {"truth": unit(r["truth"]), "low": unit(r["posterior_low"]), "median": unit(r["posterior_median"]),
                      "high": unit(r["posterior_high"]), "inside": bool(r["truth_inside"]),
                      "width_ratio": float(r["width_ratio"])}
    return pd.DataFrame(rows).T


def recovery_figure(tables: list[pd.DataFrame], titles: list[str]) -> plt.Figure:
    units = [in_prior_units(t) for t in tables]
    order = pd.concat([u["width_ratio"].astype(float) for u in units], axis=1).mean(axis=1).sort_values().index
    n = len(order)
    y = np.arange(n)[::-1]
    fig, axes = plt.subplots(1, len(units), figsize=(5.3 * len(units) + 0.6, 0.32 * n + 1.9), sharey=True, squeeze=False)
    for ax, u, title in zip(axes[0], units, titles):
        u = u.reindex(order)
        inside = u["inside"].to_numpy(dtype=bool)
        ax.axvspan(-1, 1, color=GREY_LIGHT, alpha=0.6, lw=0, zorder=0)
        ax.hlines(y, u["low"].astype(float), u["high"].astype(float), color=BLUE, lw=3, zorder=2)
        ax.scatter(u["median"].astype(float), y, s=34, color=BLUE, edgecolors="white", linewidths=1.1, zorder=3)
        truth = u["truth"].to_numpy(dtype=float)
        ax.scatter(truth[inside], y[inside], s=52, marker="D", color=ORANGE, edgecolors="white", linewidths=1.1, zorder=4)
        ax.scatter(truth[~inside], y[~inside], s=62, marker="D", facecolors="white", edgecolors=ORANGE, linewidths=1.8,
                   zorder=4)
        reach = float(np.nanmax(np.abs(u[["low", "high", "truth"]].to_numpy(dtype=float))))
        ax.set_xlim(-max(1.6, reach + 0.2), max(1.6, reach + 0.2))
        ax.set_xticks([-1, 0, 1], ["prior 5%", "median", "prior 95%"])
        ax.set_title(f"{title}: {int(inside.sum())} of {n} recovered", loc="left")
        ax.grid(False)
        ax.tick_params(axis="y", length=0)
    axes[0, 0].set_yticks(y, [talk.symbol_of(name) for name in order])
    axes[0, 0].set_ylim(-0.7, n - 0.3)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=GREY_LIGHT, alpha=0.6, label="prior 90% range"),
               plt.Line2D([], [], color=BLUE, lw=3, marker="o", mec="white", label="posterior median and 90% range"),
               plt.Line2D([], [], marker="D", ls="", color=ORANGE, mec="white", ms=8, label="true value, recovered"),
               plt.Line2D([], [], marker="D", ls="", mfc="white", mec=ORANGE, mew=1.8, ms=8, label="true value, missed")]
    fig.legend(handles=handles, loc="lower center", ncols=2 if len(units) == 1 else 4, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07 if len(units) > 1 else 0.11, 1, 1))
    return fig


# ---------------------------------------------------------------------------
# Figure AE: pairwise plot
# ---------------------------------------------------------------------------

def tied_quantities(corr: pd.DataFrame, most: int = 8, strong: float = STRONG) -> list[str]:
    """Quantities with at least one strong tie, ordered by their strongest tie, at most `most`."""
    off = corr.abs().where(~np.eye(len(corr), dtype=bool), 0.0)
    strongest = off.max(axis=1).sort_values(ascending=False)
    chosen = [name for name, value in strongest.items() if value >= strong][:most]
    if len(chosen) < 2:                       # no strong tie at all: show the most tied ones anyway
        chosen = list(strongest.index[:min(most, 4)])
    return [name for name in corr.index if name in chosen]      # keep the usual order of quantities


def pairwise_figure(draws: dict[str, np.ndarray], corr: pd.DataFrame, n_points: int = 1200) -> plt.Figure:
    names = list(draws)
    k = len(names)
    size = min(talk.FULL_WIDTH, 1.55 * k + 1.0)
    fig, axes = plt.subplots(k, k, figsize=(size, size), squeeze=False)
    keep = np.random.default_rng(0).choice(len(draws[names[0]]), size=min(len(draws[names[0]]), n_points), replace=False)
    for i, row in enumerate(names):
        for j, col in enumerate(names):
            ax = axes[i, j]
            if j > i:
                ax.axis("off")
                continue
            ax.grid(False)
            log_x = talk.on_log_axis(col)
            if i == j:                                        # the posterior of one quantity
                x, yv = talk.density(draws[col], log_x)
                if x is not None:
                    ax.fill_between(x, yv, color=BLUE, alpha=0.25, lw=0)
                    ax.plot(x, yv, color=BLUE, lw=1.6)
                ax.set_yticks([])
            else:                                             # the cloud of one pair
                r = float(corr.loc[row, col])
                if abs(r) >= STRONG:
                    ax.set_facecolor(TINT)
                ax.scatter(draws[col][keep], draws[row][keep], s=5, color=BLUE, alpha=0.35, linewidths=0)
                ax.text(0.05, 0.93, f"{r:+.2f}".replace("-", "−"), transform=ax.transAxes, va="top", fontsize=9.5,
                        color=INK_SOFT, bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.0))
                if talk.on_log_axis(row):
                    ax.set_yscale("log")
            if log_x:
                ax.set_xscale("log")
            ax.set_xticks([])                                 # the shape matters here, not the numbers
            ax.set_xticks([], minor=True)
            if i != j:
                ax.set_yticks([])
                ax.set_yticks([], minor=True)
            if i == k - 1:
                ax.set_xlabel(talk.symbol_of(col))
            if j == 0 and i > 0:
                ax.set_ylabel(talk.symbol_of(row), rotation=0, ha="right", va="center")
    handles = [plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=5, label="posterior draws"),
               plt.Rectangle((0, 0), 1, 1, facecolor=TINT, edgecolor=GREY, lw=0.8,
                             label=f"strong tie (correlation of {STRONG:g} or more)")]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.98))
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------

def recovery_table(run: str, figures: Path) -> pd.DataFrame:
    return pd.read_csv(figures / f"fake_data_recovery{run.removeprefix('fake_data_check')}.csv", index_col=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--recovery", nargs="+", default=[], help="fake-data check(s), without .nc: one figure each")
    parser.add_argument("--intervals", action="store_true",
                        help="recovery in the compact interval layout, all runs side by side in one figure")
    parser.add_argument("--pairs", default=None, help="a saved fit (real or fake data), without .nc")
    parser.add_argument("--max", type=int, default=8, help="at most this many quantities in the pairwise plot")
    args = parser.parse_args()
    config = load_config()
    apply_style()
    plt.rcParams.update(talk.TALK_STYLE)
    figures = resolve_path(config["paths"]["figures"])

    if args.recovery and args.intervals:
        tables = [recovery_table(run, figures) for run in args.recovery]
        titles = [f"true forest {2 if '_truth2' in run else 1}" for run in args.recovery]
        save_figure(recovery_figure(tables, titles), f"fig_AD_{args.recovery[0]}_recovery_intervals", config)
    for run in ([] if args.intervals else args.recovery):
        idata = talk.load(run, config)
        names = [n for n in talk.parameter_names(idata.posterior) if n != "sigma_m"]
        posterior = {n: talk.chains_of(idata.posterior, n).reshape(-1) for n in names}
        fig = recovery_density_figure(talk.prior_draws(idata, names, config), posterior, recovery_table(run, figures))
        save_figure(fig, f"fig_AD_{run}_recovery", config)

    if args.pairs:
        post = talk.load(args.pairs, config).posterior
        names = [n for n in talk.parameter_names(post) if n != "sigma_m"]
        corr = talk.correlation_matrix(post, names)
        chosen = tied_quantities(corr, most=args.max)
        print("  pairwise plot of: " + ", ".join(chosen))
        draws = {name: talk.chains_of(post, name).reshape(-1) for name in chosen}
        save_figure(pairwise_figure(draws, corr.loc[chosen, chosen]), f"fig_AE_{args.pairs}_pairwise", config)


if __name__ == "__main__":
    main()
