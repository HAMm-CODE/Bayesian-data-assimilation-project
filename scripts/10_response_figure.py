"""Which quantities respond to the data: prior against posterior for every quantity
of one fake-data check, sorted by how much NEE narrowed each one (for talks; D63).

The same panels as figure T of talk_figures.py (grey = prior, blue = posterior,
orange = the true value), but ordered by the width ratio of the recovery table
(posterior 90% width divided by prior 90% width; 1 means NEE taught nothing).
Panels at or above the cut (D49: 0.5) get a light purple background: these are
the quantities that may be held fixed. Each panel shows its width ratio; with
--also, the ratio from a second true forest is shown next to it, so one can see
whether the pattern holds for both truths. Nothing is sampled.

Run from the repository root, for example:
    python scripts/10_response_figure.py fake_data_check_cold_dense_ta995_d2000 \
        --also fake_data_check_cold_truth2_dense_ta995_d2000

Saves (in figures/, pdf and png):
    fig_T2_<run>_by_response
and prints the ordered table.
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path
from dalec.plotting import BLUE, GREY, GREY_LIGHT, ORANGE, apply_style, save_figure

_spec = importlib.util.spec_from_file_location(
    "talk_figures", Path(__file__).resolve().parent / "talk_figures.py")
talk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(talk)

TINT = "#F1EAF8"      # light purple: may be fixed
EDGE = "#4E008E"      # Tampere purple


def recovery_table(run: str, config) -> pd.DataFrame:
    path = resolve_path(config["paths"]["figures"]) / f"fake_data_recovery{run.removeprefix('fake_data_check')}.csv"
    return pd.read_csv(path, index_col=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run", help="fake-data check in results/, without .nc")
    parser.add_argument("--also", default=None, help="the same model's check with a second true forest")
    parser.add_argument("--cut", type=float, default=0.5, help="width ratio at or above which a quantity may be fixed (D49)")
    args = parser.parse_args()

    config = load_config()
    apply_style()
    plt.rcParams.update(talk.TALK_STYLE)
    idata = talk.load(args.run, config)
    post = idata.posterior
    names = [n for n in talk.parameter_names(post) if n != "sigma_m"]
    table = recovery_table(args.run, config)
    ratio = table["width_ratio"].astype(float)
    ratio2 = recovery_table(args.also, config)["width_ratio"].astype(float) if args.also else None
    truth = {n: float(v) for n, v in table["truth"].items()}

    key = ratio if ratio2 is None else (ratio + ratio2.reindex(ratio.index)) / 2
    names = sorted(names, key=lambda n: -float(key.get(n, -1.0)))
    prior = talk.prior_draws(idata, names, config)

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
        x, y = talk.density(talk.chains_of(post, name).reshape(-1), log_axis)
        if x is not None:
            ax.plot(x, y, color=BLUE, lw=2)
        if name in truth:
            ax.axvline(truth[name], color=ORANGE, lw=2)
        if log_axis:
            ax.set_xscale("log")
        exponent = talk.scale_exponent(talk.chains_of(post, name)) if log_axis else None
        talk.few_ticks(ax, "x", 2, log_axis, exponent)
        ax.set_yticks([])
        ax.set_ylim(0, 1.32)
        ax.set_title(talk.with_scale(talk.short_label(name), exponent))
        text = f"{ratio[name]:.2f}" if name in ratio else ""
        if ratio2 is not None and name in ratio2:
            text += f" / {ratio2[name]:.2f}"
        ax.text(0.97, 0.96, text, transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color="#2b2b33",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.2), zorder=5)
        if float(key.get(name, 0.0)) >= args.cut:
            ax.set_facecolor(TINT)
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color(EDGE)
                spine.set_linewidth(1.2)

    handles = [plt.Line2D([], [], color=GREY, lw=7, alpha=0.4, label="prior"),
               plt.Line2D([], [], color=BLUE, lw=2.5, label="posterior"),
               plt.Line2D([], [], color=ORANGE, lw=2.5, label="true value (fake data)"),
               plt.Rectangle((0, 0), 1, 1, facecolor=TINT, edgecolor=EDGE, lw=1.2,
                             label=f"width ratio {args.cut:g} or more: may be fixed")]
    fig.legend(handles=handles, loc="lower center", ncols=len(handles), bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_figure(fig, f"fig_T2_{args.run}_by_response", config)

    out = pd.DataFrame({"width_ratio": ratio.reindex(names)})
    if ratio2 is not None:
        out["width_ratio_truth2"] = ratio2.reindex(names)
    out["may_be_fixed"] = key.reindex(names) >= args.cut
    print(out.round(3).to_string())


if __name__ == "__main__":
    main()
