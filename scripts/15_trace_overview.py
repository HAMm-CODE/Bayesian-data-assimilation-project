"""Trace plots of every sampled quantity, for the seminar. Two outputs.

  1. fig_R1_<fit>_trace_overview: every quantity on ONE page, density and trace side by side (small panels).
  2. fig_R2_<fit>_trace_density_all.pdf: ONE PDF with every quantity, density per
     chain next to its trace, 6 quantities per page (the layout of figure R in
     talk_figures.py, but all pages in a single file and with the notation-sheet labels).


Nothing is sampled; the script only reads the saved posterior. Each small panel
shows the draws of one quantity, in order, one colour per chain. If the sampler
worked, the four chains overlap and look like one fuzzy band ("hairy
caterpillar") that stays at one level. The densities are not repeated here: the
prior and posterior densities are in figure T (talk_figures.py).

Labels follow the notation sheet (D70): f_fol, c_lma, c_ronset, C^0_lab, and so on.

Run from the repository root, for example:
    python scripts/15_trace_overview.py fit_model_error_MC_cold_dense_ta995_d2000

Saves in figures/:
    fig_R1_<fit>_trace_overview        (pdf and png)
    fig_R2_<fit>_trace_density_all.pdf (all pages in one PDF)
"""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

from dalec.config import load_config, resolve_path
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
COLUMNS = 4


def thesis_label(name: str) -> str:
    return THESIS_LABELS.get(name, talk.short_label(name))


def trace_overview_figure(post, names: list[str], columns: int = COLUMNS) -> plt.Figure:
    """Every quantity on ONE page: a pair of panels each, the density of each chain (narrow)
    next to its trace (wide). `columns` pairs per row."""
    n_chains = post[names[0]].shape[0]
    rows = int(np.ceil(len(names) / columns))
    fig, axes = plt.subplots(rows, 2 * columns, figsize=(talk.FULL_WIDTH, 1.15 * rows + 0.8), squeeze=False,
                             gridspec_kw={"width_ratios": [1, 1.7] * columns})
    for slot in range(len(names), rows * columns):
        axes[slot // columns, 2 * (slot % columns)].axis("off")
        axes[slot // columns, 2 * (slot % columns) + 1].axis("off")
    for slot, name in enumerate(names):
        row, col = divmod(slot, columns)
        dens, trace = axes[row, 2 * col], axes[row, 2 * col + 1]
        values = talk.chains_of(post, name)
        log_axis = talk.on_log_axis(name)
        for c in range(n_chains):
            colour = talk.CHAIN_COLOURS[c % len(talk.CHAIN_COLOURS)]
            x, y = talk.density(values[c], log_axis)
            if x is not None:
                dens.plot(x, y, color=colour, lw=1.4)
            trace.plot(values[c], color=colour, lw=0.45, alpha=0.75)
        if log_axis:
            dens.set_xscale("log")
            trace.set_yscale("log")
        exponent = talk.scale_exponent(values) if log_axis else None
        talk.few_ticks(dens, "x", 2, log_axis, exponent)
        talk.few_ticks(trace, "y", 3, log_axis, exponent)
        dens.tick_params(axis="x", labelsize=8)
        dens.set_yticks([])
        dens.set_title(talk.with_scale(thesis_label(name), exponent), loc="left")
        trace.set_xlim(0, values.shape[1])
        trace.set_xticks([0, values.shape[1]])
        dens.grid(False)
        trace.grid(False)
        last_in_column = slot + columns >= len(names)
        if not last_in_column:
            trace.tick_params(labelbottom=False)
        else:
            trace.set_xlabel("draw")
    talk.bottom_legend(fig, n_chains)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig


def density_trace_pages(post, names: list[str], per_page: int = 6) -> list[plt.Figure]:
    """Pages of two quantities per row: the density of each chain (narrow panel), then its draws in order."""
    n_chains = post[names[0]].shape[0]
    figures = []
    for start in range(0, len(names), per_page):
        page = names[start:start + per_page]
        rows = int(np.ceil(len(page) / 2))
        fig, axes = plt.subplots(rows, 4, figsize=(talk.FULL_WIDTH, 1.6 * rows + 0.8),
                                 gridspec_kw={"width_ratios": [1, 1.9, 1, 1.9]}, squeeze=False)
        for slot in range(len(page), 2 * rows):                  # an odd last page
            axes[slot // 2, 2 * (slot % 2)].axis("off")
            axes[slot // 2, 2 * (slot % 2) + 1].axis("off")
        for slot, name in enumerate(page):
            dens, trace = axes[slot // 2, 2 * (slot % 2)], axes[slot // 2, 2 * (slot % 2) + 1]
            values = talk.chains_of(post, name)
            log_axis = talk.on_log_axis(name)
            for c in range(n_chains):
                colour = talk.CHAIN_COLOURS[c % len(talk.CHAIN_COLOURS)]
                x, y = talk.density(values[c], log_axis)
                if x is not None:
                    dens.plot(x, y, color=colour, lw=1.6)
                trace.plot(values[c], color=colour, lw=0.6, alpha=0.75)
            if log_axis:
                dens.set_xscale("log")
                trace.set_yscale("log")
            exponent = talk.scale_exponent(values) if log_axis else None
            talk.few_ticks(dens, "x", 2, log_axis, exponent)
            talk.few_ticks(trace, "y", 3, log_axis, exponent)
            dens.set_yticks([])
            dens.set_title(talk.with_scale(thesis_label(name), exponent), loc="left")
            trace.set_xlim(0, values.shape[1])
            trace.xaxis.set_major_locator(MaxNLocator(4, integer=True))
            if slot + 2 < len(page):                             # draw numbers on the last row only
                trace.tick_params(labelbottom=False)
            else:
                trace.set_xlabel("draw")
        talk.bottom_legend(fig, n_chains)
        fig.tight_layout(rect=(0, 0.07, 1, 1))
        figures.append(fig)
    return figures


def save_one_pdf(figures: list[plt.Figure], path: Path) -> None:
    """All pages in a single PDF file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(path) as pdf:
        for fig in figures:
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
    print(f"  saved {path} ({len(figures)} pages)")


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
    out_dir = resolve_path(config["paths"]["figures"])
    save_one_pdf(density_trace_pages(post, names), out_dir / f"fig_R2_{args.fit}_trace_density_all.pdf")


if __name__ == "__main__":
    main()
