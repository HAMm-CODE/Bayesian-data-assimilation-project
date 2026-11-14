"""Trace plots of every sampled quantity on ONE slide (for the seminar).

Nothing is sampled; the script only reads the saved posterior. Each small panel
shows the draws of one quantity, in order, one colour per chain. If the sampler
worked, the four chains overlap and look like one fuzzy band ("hairy
caterpillar") that stays at one level. The densities are not repeated here: the
prior and posterior densities are in figure T (talk_figures.py).

Labels follow the notation sheet (D70): f_fol, c_lma, c_ronset, C^0_lab, and so on.

Run from the repository root, for example:
    python scripts/15_trace_overview.py fit_model_error_MC_cold_dense_ta995_d2000

Saves (in figures/, pdf and png):
    fig_R1_<fit>_trace_overview
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from dalec.config import load_config
from dalec.plotting import apply_style, save_figure

_spec = importlib.util.spec_from_file_location("talk_figures", Path(__file__).resolve().parent / "talk_figures.py")
talk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(talk)

# Names as in the notation sheet. The code names stay in the code.
THESIS_LABELS = {
    "f_auto": r"$f_\mathrm{auto}$",
    "frac_fol": r"$f_\mathrm{fol}$", "frac_lab": r"$f_\mathrm{lab}$", "frac_roo": r"$f_\mathrm{roo}$",
    "theta_roo": r"$\theta_\mathrm{roo}$", "theta_woo": r"$\theta_\mathrm{woo}$",
    "theta_lit": r"$\theta_\mathrm{lit}$", "theta_som": r"$\theta_\mathrm{som}$",
    "theta_min": r"$\theta_\mathrm{min}$", "temperature_exponent": r"$\Theta$",
    "d_onset": r"$d_\mathrm{onset}$", "cr_onset": r"$c_\mathrm{ronset}$",
    "d_fall": r"$d_\mathrm{fall}$", "cr_fall": r"$c_\mathrm{rfall}$",
    "c_lf": r"$c_\mathrm{lf}$", "lma": r"$c_\mathrm{lma}$", "ceff": r"$c_\mathrm{eff}$",
    "c_lab_0": r"$C^0_\mathrm{lab}$", "c_fol_0": r"$C^0_\mathrm{fol}$", "c_roo_0": r"$C^0_\mathrm{roo}$",
    "c_woo_0": r"$C^0_\mathrm{woo}$", "c_lit_0": r"$C^0_\mathrm{lit}$", "c_som_0": r"$C^0_\mathrm{som}$",
    "cold_t_shutdown": r"$T_\mathrm{shutdown}$", "cold_t_full": r"$T_\mathrm{full}$",
    "sigma_m": r"$\sigma_m$",
}
COLUMNS = 6


def thesis_label(name: str) -> str:
    return THESIS_LABELS.get(name, talk.short_label(name))


def trace_overview_figure(post, names: list[str], columns: int = COLUMNS) -> plt.Figure:
    """One small trace panel per quantity, `columns` panels per row."""
    n_chains = post[names[0]].shape[0]
    rows = int(np.ceil(len(names) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(talk.FULL_WIDTH, 1.35 * rows + 0.8), squeeze=False)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    last_row_of = {c: max(r for r in range(rows) if r * columns + c < len(names)) for c in range(columns)
                   if c < len(names)}
    for i, (ax, name) in enumerate(zip(axes.flat, names)):
        values = talk.chains_of(post, name)
        log_axis = talk.on_log_axis(name)
        for c in range(n_chains):
            ax.plot(values[c], color=talk.CHAIN_COLOURS[c % len(talk.CHAIN_COLOURS)], lw=0.5, alpha=0.75)
        if log_axis:
            ax.set_yscale("log")
        exponent = talk.scale_exponent(values) if log_axis else None
        talk.few_ticks(ax, "y", 3, log_axis, exponent)
        ax.set_xlim(0, values.shape[1])
        ax.set_xticks([0, values.shape[1] // 2, values.shape[1]])
        ax.set_title(talk.with_scale(thesis_label(name), exponent), loc="left")
        ax.grid(False)
        row, col = divmod(i, columns)
        if row != last_row_of[col]:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("draw")
    talk.bottom_legend(fig, n_chains)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fit", help="real-data fit, without .nc")
    args = parser.parse_args()
    config = load_config()
    apply_style()
    plt.rcParams.update(talk.TALK_STYLE)
    post = talk.load(args.fit, config).posterior
    names = talk.parameter_names(post)
    print(f"  {len(names)} quantities: {', '.join(names)}")
    save_figure(trace_overview_figure(post, names), f"fig_R1_{args.fit}_trace_overview", config)


if __name__ == "__main__":
    main()
