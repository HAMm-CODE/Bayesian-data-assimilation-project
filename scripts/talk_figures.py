"""Figures and tables for talks (and the thesis) from one saved posterior.

Works on any fit saved by scripts 04 or 05 (results/<label>.nc). Nothing is
sampled; the script only reads the posterior and, with --predict, runs the
numpy model forward for posterior draws.

Run from the repository root, for example:
    python scripts/talk_figures.py fake_data_check_dense_ta95 --truth figures/fake_data_recovery_dense_ta95.csv
    python scripts/talk_figures.py fit_model_error_cold_quick --predict

Makes (in figures/, pdf and png), sized for a 16:9 slide:
    R  fig_R_<label>_trace_p<k>          every parameter: density per chain and the draws in
                                         order (the "hairy caterpillar"), 6 per page
    S  fig_S_<label>_rank_p<k>           every parameter: rank plot per chain, 12 per page
    T  fig_T_<label>_prior_posterior     every parameter: prior (grey) against posterior (blue),
                                         true value in orange for fake-data runs (RQ1)
    U  fig_U_<label>_correlations        posterior correlation of every pair (RQ3)
       fig_U_<label>_pairs               the six most strongly correlated pairs (RQ3)
  with --predict (real-data fits only):
    V  fig_V_<label>_predictive_nee      50 posterior runs of daily NEE over the observed
                                         ("spaghetti" plot), first two years of the block
    W  fig_W_<label>_pools               the six hidden pools at the start of each year,
                                         median and 90% band (RQ1, RQ2)
    X  fig_X_<label>_uncertainty         where the uncertainty comes from, month by month:
                                         parameters, measurement error, model error (RQ2)
Tables (in figures/):
    <label>_summary.csv                  mean, sd, 5% and 95%, R-hat and ESS for every parameter
    <label>_summary_part1.tex, _part2.tex   the same as two booktabs tables for slides

The figures carry no titles: on a slide the frame title says what they show,
and in the thesis the caption does.
"""

import argparse
from typing import Any

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FixedFormatter, FixedLocator, MaxNLocator, NullLocator
from scipy import stats

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.model_numpy import POOLS, drivers_from_frame, run
from dalec.plotting import (
    BLUE, GREY, GREY_LIGHT, INK, INK_SOFT, MONTHS, ORANGE, UNIT_FLUX, apply_style, heading,
    save_figure,
)
from dalec.posterior import QUANTITIES, posterior_parameters
from dalec.priors import (
    ACM_PRIORS, COLD_PRIORS, PRIORS, priors_from_attrs, sample_acm, sample_cold, sample_prior,
    to_prior_scale,
)

# One colour per chain. Checked as a set for colour-blind safety, every pair
# compared (dataviz validator, light mode, --pairs all); a legend names them.
CHAIN_COLOURS = ("#2a78d6", "#eda100", "#e87ba4", "#008300")

# Figures are drawn about twice as large as they appear on a 16:9 slide, with
# large fonts, so the text comes out at 5 to 6 pt on the slide (the size of the
# deck's small text). FULL_WIDTH fills a slide, HALF_WIDTH one of two columns.
FULL_WIDTH = 12.0
HALF_WIDTH = 6.6
TALK_STYLE = {"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 12,
              "xtick.labelsize": 10.5, "ytick.labelsize": 10.5, "legend.fontsize": 12}
RANK_PER_PAGE = 12

# Plot labels (matplotlib mathtext, which LaTeX also reads) and units.
LABELS = {
    "f_auto": r"$f_\mathrm{auto}$", "f_lab": r"$f_\mathrm{lab}$", "f_fol": r"$f_\mathrm{fol}$",
    "f_roo": r"$f_\mathrm{roo}$", "f_woo": r"$f_\mathrm{woo}$",
    "frac_fol": r"$s_\mathrm{fol}$", "frac_lab": r"$s_\mathrm{lab}$", "frac_roo": r"$s_\mathrm{roo}$",
    "theta_roo": r"$\theta_\mathrm{roo}$ (d$^{-1}$)", "theta_woo": r"$\theta_\mathrm{woo}$ (d$^{-1}$)",
    "theta_lit": r"$\theta_\mathrm{lit}$ (d$^{-1}$)", "theta_som": r"$\theta_\mathrm{som}$ (d$^{-1}$)",
    "theta_min": r"$\theta_{\mathrm{min}}$ (d$^{-1}$)", "temperature_exponent": r"$\Theta$ (K$^{-1}$)",
    "d_onset": r"$d_\mathrm{onset}$ (day)", "cr_onset": r"$cr_\mathrm{onset}$ (d)",
    "d_fall": r"$d_\mathrm{fall}$ (day)", "cr_fall": r"$cr_\mathrm{fall}$ (d)",
    "c_lf": r"$c_\mathrm{lf}$ (yr$^{-1}$)", "lma": r"LMA (g C m$^{-2}$)", "ceff": r"$c_\mathrm{eff}$",
    "c_lab_0": r"$C_\mathrm{lab,0}$", "c_fol_0": r"$C_\mathrm{fol,0}$", "c_roo_0": r"$C_\mathrm{roo,0}$",
    "c_woo_0": r"$C_\mathrm{woo,0}$", "c_lit_0": r"$C_\mathrm{lit,0}$", "c_som_0": r"$C_\mathrm{som,0}$",
    "cold_t_shutdown": r"$T_\mathrm{shutdown}$ (°C)", "cold_t_full": r"$T_\mathrm{full}$ (°C)",
    "acm_a7": r"$a_7$ (ACM light-use efficiency)", "sigma_m": r"$\sigma_m$ (" + UNIT_FLUX + ")",
}
# Where the long label does not fit (small panels)
SHORT_LABELS = {"theta_min": r"$\theta_{\mathrm{min}}$", "acm_a7": r"$a_7$", "sigma_m": r"$\sigma_m$"}
POOL_LABELS = {"c_lab": "labile", "c_fol": "foliage", "c_roo": "roots",
               "c_woo": "wood", "c_lit": "litter", "c_som": "soil organic matter"}


def label_of(name: str) -> str:
    return LABELS.get(name, name)


def short_label(name: str) -> str:
    return SHORT_LABELS.get(name, label_of(name))


def symbol_of(name: str) -> str:
    """The symbol alone, without units (tables and the correlation axes)."""
    return label_of(name).split(" (")[0]


def on_log_axis(name: str) -> bool:
    """Rates, pools and other log-normal quantities are easier to read on a log axis."""
    return (name in PRIORS and PRIORS[name][0] == "lognormal") or name in ACM_PRIORS


# ---------------------------------------------------------------------------
# Axis ticks that stay readable in small panels
# ---------------------------------------------------------------------------

def scale_exponent(values: np.ndarray) -> int | None:
    """Power of ten taken out of the tick labels of a very small quantity, so a
    rate of 0.00003 per day reads "3" on an axis labelled (10^-5 d^-1)."""
    middle = float(np.median(np.abs(np.asarray(values, dtype=float))))
    if middle == 0 or middle >= 1e-2:
        return None
    return int(np.floor(np.log10(middle)))


def with_scale(label: str, exponent: int | None) -> str:
    """Add the power of ten to a label: 'θ (d⁻¹)' becomes 'θ (10⁻⁵ d⁻¹)'."""
    if exponent is None:
        return label
    power = rf"10$^{{{exponent}}}$"
    if " (" in label:
        symbol, unit = label.split(" (", 1)
        return f"{symbol} ({power} {unit}"
    return f"{label} ({power})"


def tick_text(value: float) -> str:
    """Plain numbers (0.05, 300, 4000); powers of ten only for very small or large ones."""
    if value == 0:
        return "0"
    if 1e-3 <= abs(value) < 1e5:
        return f"{value:.2g}" if abs(value) < 1 else f"{value:g}"
    mantissa, exponent = f"{value:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def log_ticks(lo: float, hi: float, n: int) -> list[float]:
    """Up to n ticks between lo and hi, away from the edges: 1, 2 and 5 times a
    power of ten when enough fall inside, otherwise evenly spread values rounded
    to two figures. Spread out so their labels do not run into each other."""
    margin = (hi / lo) ** 0.06
    inner_lo, inner_hi = lo * margin, hi / margin
    nice = [m * 10.0 ** k for k in range(int(np.floor(np.log10(lo))) - 1, int(np.ceil(np.log10(hi))) + 1)
            for m in (1, 2, 5)]
    inside = [t for t in nice if inner_lo <= t <= inner_hi]
    if len(inside) >= min(n, 2):                 # two round numbers beat three odd ones
        picks = sorted(set(np.round(np.linspace(0, len(inside) - 1, n)).astype(int)))
        return [inside[i] for i in picks]
    positions = np.exp(np.linspace(np.log(lo), np.log(hi), n + 2)[1:-1])
    return sorted({float(format(p, ".2g")) for p in positions})


def tidy_log_axis(ax, which: str, n: int = 3, exponent: int | None = None) -> None:
    """At most n plain-number ticks on a log axis. Matplotlib's default labels
    every minor tick of a short log axis, and the labels overlap. With an
    exponent the labels are divided by 10**exponent (see with_scale)."""
    lo, hi = sorted(ax.get_xlim() if which == "x" else ax.get_ylim())
    axis = ax.xaxis if which == "x" else ax.yaxis
    ticks = log_ticks(lo, hi, n)
    scale = 1.0 if exponent is None else 10.0 ** exponent
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(FixedFormatter([tick_text(float(f"{t / scale:.3g}")) for t in ticks]))
    axis.set_minor_locator(NullLocator())


def few_ticks(ax, which: str, n: int, log_axis: bool, exponent: int | None = None) -> None:
    """About n ticks on one axis, log or linear."""
    if log_axis:
        tidy_log_axis(ax, which, n, exponent)
    else:
        (ax.xaxis if which == "x" else ax.yaxis).set_major_locator(MaxNLocator(n))


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load(label: str, config: dict[str, Any]):
    path = resolve_path(config["paths"]["results"]) / f"{label}.nc"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found; the label is the file name without .nc")
    return az.from_netcdf(path)


def parameter_names(post) -> list[str]:
    """Every sampled quantity in the fit, in the D34 to D39 order, then the extras.
    Quantities held fixed (D54) are left out: a constant has no distribution."""
    candidates = list(QUANTITIES) + list(COLD_PRIORS) + list(ACM_PRIORS) + ["sigma_m"]
    names = []
    for name in candidates:
        if name in post.data_vars and float(np.ptp(np.asarray(post[name].values))) > 0:
            names.append(name)
    return names


def chains_of(post, name: str) -> np.ndarray:
    """(chains, draws) array of one quantity."""
    return np.asarray(post[name].values, dtype=float).reshape(post[name].shape[0], -1)


def prior_draws(idata, names: list[str], config: dict[str, Any], n: int = 20_000) -> dict[str, np.ndarray]:
    """Draws from the priors the fit was made with (saved in the file since D57)."""
    rng = np.random.default_rng(1)
    _, raw = sample_prior(n, rng, priors_from_attrs(idata.attrs))
    raw.update(sample_cold(n, rng))
    raw.update(sample_acm(n, rng))
    raw["sigma_m"] = np.abs(rng.standard_normal(n)) * float(config["likelihood"]["model_error_prior_sd"])
    return {name: raw[name] for name in names if name in raw}


def read_truth(path: str | None) -> dict[str, float]:
    if not path:
        return {}
    table = pd.read_csv(path, index_col=0)
    return {name: float(value) for name, value in table["truth"].items()}


# ---------------------------------------------------------------------------
# Densities on the scale each quantity is read on
# ---------------------------------------------------------------------------

def density(values: np.ndarray, log_axis: bool, grid_size: int = 200):
    """Kernel density curve, scaled to a peak of 1. On a log axis the density
    is computed for log(value), which is what the eye reads on that axis."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    x = np.log(values) if log_axis else values
    if np.ptp(x) == 0:
        return None, None
    kde = stats.gaussian_kde(x)
    lo, hi = np.percentile(x, [0.5, 99.5])
    pad = 0.1 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, grid_size)
    y = kde(grid)
    return (np.exp(grid) if log_axis else grid), y / y.max()


# ---------------------------------------------------------------------------
# R and S: every parameter, chain by chain
# ---------------------------------------------------------------------------

def pages(names: list[str], per_page: int) -> list[list[str]]:
    return [names[i:i + per_page] for i in range(0, len(names), per_page)]


def bottom_legend(fig, n_chains: int, extra: tuple = ()) -> None:
    """One entry per chain, plus any reference lines, in a row under the panels."""
    handles = [plt.Line2D([], [], color=CHAIN_COLOURS[c % len(CHAIN_COLOURS)], lw=2.5,
                          label=f"chain {c + 1}") for c in range(n_chains)]
    handles += list(extra)
    fig.legend(handles=handles, loc="lower center", ncols=len(handles), bbox_to_anchor=(0.5, 0.0))


def trace_figures(post, names, label, config, truth, per_page) -> None:
    """Two parameters per row: the density of each chain (narrow panel), then its draws in order."""
    n_chains = post[names[0]].shape[0]
    extra = (plt.Line2D([], [], color=ORANGE, lw=2, ls="--", label="true value"),) if truth else ()
    for k, page in enumerate(pages(names, per_page), start=1):
        rows = int(np.ceil(len(page) / 2))
        fig, axes = plt.subplots(rows, 4, figsize=(FULL_WIDTH, 1.6 * rows + 0.8),
                                 gridspec_kw={"width_ratios": [1, 1.9, 1, 1.9]}, squeeze=False)
        for slot in range(len(page), 2 * rows):                  # an odd last page
            axes[slot // 2, 2 * (slot % 2)].axis("off")
            axes[slot // 2, 2 * (slot % 2) + 1].axis("off")
        for slot, name in enumerate(page):
            dens, trace = axes[slot // 2, 2 * (slot % 2)], axes[slot // 2, 2 * (slot % 2) + 1]
            values = chains_of(post, name)
            log_axis = on_log_axis(name)
            for c in range(n_chains):
                colour = CHAIN_COLOURS[c % len(CHAIN_COLOURS)]
                x, y = density(values[c], log_axis)
                if x is not None:
                    dens.plot(x, y, color=colour, lw=1.6)
                trace.plot(values[c], color=colour, lw=0.6, alpha=0.75)
            if name in truth:
                dens.axvline(truth[name], color=ORANGE, lw=1.8, ls="--")
                trace.axhline(truth[name], color=ORANGE, lw=1.8, ls="--")
            if log_axis:
                dens.set_xscale("log")
                trace.set_yscale("log")
            exponent = scale_exponent(values) if log_axis else None
            few_ticks(dens, "x", 2, log_axis, exponent)
            few_ticks(trace, "y", 3, log_axis, exponent)
            dens.set_yticks([])
            dens.set_title(with_scale(short_label(name), exponent), loc="left")
            trace.set_xlim(0, values.shape[1])
            trace.xaxis.set_major_locator(MaxNLocator(4, integer=True))
            if slot + 2 < len(page):                             # draw numbers on the last row only
                trace.tick_params(labelbottom=False)
            else:
                trace.set_xlabel("draw")
        bottom_legend(fig, n_chains, extra)
        fig.tight_layout(rect=(0, 0.07, 1, 1))
        save_figure(fig, f"fig_R_{label}_trace_p{k}", config)


def rank_figures(post, names, label, config, per_page: int = RANK_PER_PAGE, bins: int = 20) -> None:
    """Rank plots (Vehtari et al. 2021): all draws ranked together; if the chains
    agree, every chain's ranks are spread evenly, so its bars sit flat."""
    n_chains, n_draws = post[names[0]].shape[0], post[names[0]].shape[1]
    expected = n_draws / bins
    extra = (plt.Line2D([], [], color=INK_SOFT, lw=1.2, ls="--", label="flat if the chains agree"),)
    for k, page in enumerate(pages(names, per_page), start=1):
        rows = int(np.ceil(len(page) / 3))
        fig, axes = plt.subplots(rows, 3, figsize=(FULL_WIDTH, 1.25 * rows + 0.8), squeeze=False)
        for ax in axes.flat[len(page):]:
            ax.axis("off")
        for ax, name in zip(axes.flat, page):
            values = chains_of(post, name)
            ranks = stats.rankdata(values.reshape(-1)).reshape(values.shape)
            edges = np.linspace(0.5, values.size + 0.5, bins + 1)
            width = 1.0 / n_chains
            for c in range(n_chains):
                counts, _ = np.histogram(ranks[c], bins=edges)
                ax.bar(np.arange(bins) + c * width, counts, width=width * 0.85, align="edge",
                       color=CHAIN_COLOURS[c % len(CHAIN_COLOURS)], linewidth=0)
            ax.axhline(expected, color=INK_SOFT, lw=1.2, ls="--")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            ax.set_title(short_label(name))
        bottom_legend(fig, n_chains, extra)
        fig.tight_layout(rect=(0, 0.08, 1, 1))
        save_figure(fig, f"fig_S_{label}_rank_p{k}", config)


# ---------------------------------------------------------------------------
# T: prior against posterior, every parameter (RQ1)
# ---------------------------------------------------------------------------

def prior_posterior_figure(post, names, prior, label, config, truth) -> None:
    cols = 7 if len(names) > 20 else 5
    rows = int(np.ceil(len(names) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(FULL_WIDTH, 1.25 * rows + 0.75), squeeze=False)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    for ax, name in zip(axes.flat, names):
        log_axis = on_log_axis(name)
        if name in prior:
            x, y = density(prior[name], log_axis)
            if x is not None:
                ax.fill_between(x, y, color=GREY_LIGHT, lw=0)
                ax.plot(x, y, color=GREY, lw=1)
        x, y = density(chains_of(post, name).reshape(-1), log_axis)
        if x is not None:
            ax.plot(x, y, color=BLUE, lw=2)
        if name in truth:
            ax.axvline(truth[name], color=ORANGE, lw=2)
        if log_axis:
            ax.set_xscale("log")
        exponent = scale_exponent(chains_of(post, name)) if log_axis else None
        few_ticks(ax, "x", 2, log_axis, exponent)
        ax.set_yticks([])
        ax.set_title(with_scale(short_label(name), exponent))
    handles = [plt.Line2D([], [], color=GREY, lw=7, alpha=0.4, label="prior"),
               plt.Line2D([], [], color=BLUE, lw=2.5, label="posterior")]
    if truth:
        handles.append(plt.Line2D([], [], color=ORANGE, lw=2.5, label="true value (fake data)"))
    fig.legend(handles=handles, loc="lower center", ncols=len(handles), bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_figure(fig, f"fig_T_{label}_prior_posterior", config)


# ---------------------------------------------------------------------------
# U: which parameters move together (RQ3)
# ---------------------------------------------------------------------------

def correlation_matrix(post, names) -> pd.DataFrame:
    """Posterior correlations on each quantity's prior scale (log or logit, D32)."""
    columns = {}
    for name in names:
        values = chains_of(post, name).reshape(-1)
        if name in COLD_PRIORS or name in ACM_PRIORS or name in PRIORS or name.startswith("f_"):
            values = to_prior_scale(name, values)
        columns[name] = values
    return pd.DataFrame(columns).corr()


def correlation_figures(post, names, label, config, n_pairs: int = 6) -> None:
    corr = correlation_matrix(post, names)
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list(
        "diverging", ["#1c5cab", "#f0efec", "#e34948"])
    fig, ax = plt.subplots(figsize=(HALF_WIDTH + 0.4, 6.2))
    image = ax.imshow(corr.values, cmap=cmap, vmin=-1, vmax=1)
    symbols = [symbol_of(n) for n in names]
    ax.set_xticks(range(len(names)), symbols, rotation=90)
    ax.set_yticks(range(len(names)), symbols)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    bar = fig.colorbar(image, ax=ax, shrink=0.75, pad=0.02, ticks=[-1, -0.5, 0, 0.5, 1])
    bar.set_label("correlation after the fit")
    fig.tight_layout()
    save_figure(fig, f"fig_U_{label}_correlations", config)

    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pairs.append((abs(corr.iloc[i, j]), names[i], names[j], corr.iloc[i, j]))
    pairs = sorted(pairs, reverse=True)[:n_pairs]
    fig, axes = plt.subplots(2, 3, figsize=(FULL_WIDTH, 5.8), squeeze=False)
    for ax, (_, a, b, r) in zip(axes.flat, pairs):
        x = chains_of(post, a).reshape(-1)
        y = chains_of(post, b).reshape(-1)
        keep = np.random.default_rng(0).choice(x.size, size=min(x.size, 1500), replace=False)
        ax.scatter(x[keep], y[keep], s=14, color=BLUE, alpha=0.35, linewidths=0)
        if on_log_axis(a):
            ax.set_xscale("log")
        if on_log_axis(b):
            ax.set_yscale("log")
        exponent_a = scale_exponent(x) if on_log_axis(a) else None
        exponent_b = scale_exponent(y) if on_log_axis(b) else None
        few_ticks(ax, "x", 3, on_log_axis(a), exponent_a)
        few_ticks(ax, "y", 3, on_log_axis(b), exponent_b)
        ax.set_xlabel(with_scale(short_label(a), exponent_a))
        ax.set_ylabel(with_scale(short_label(b), exponent_b))
        ax.set_title(f"correlation {r:+.2f}".replace("-", "−"))
    fig.tight_layout()
    save_figure(fig, f"fig_U_{label}_pairs", config)


# ---------------------------------------------------------------------------
# V, W, X: running the posterior forward (real-data fits)
# ---------------------------------------------------------------------------

def fit_block(label: str, config: dict[str, Any]) -> pd.DataFrame:
    """The days the fit used: the calibration block, or its first years for quick fits."""
    data = load_model_block(config, block="calibration")
    if label.endswith("_quick"):
        years = data.index.year
        data = data[years < years[0] + int(config["sampler"]["quick"]["years"])]
    return data


def forward(idata, data, config, label, n_draws: int, seed: int):
    """Run the numpy model for n posterior draws. Returns the model output
    (daily NEE, pools at the start of each year, ...) and sigma_m per draw."""
    drivers = drivers_from_frame(data, config)
    frost_cfg = config["model"]["frost_cutoff"]
    frost_on = "_frost" in label or bool(frost_cfg["enabled"])
    frost = float(frost_cfg["threshold_degc"]) if frost_on else None
    params, sigma_m, cold, a7 = posterior_parameters(idata, n_draws, seed)
    with np.errstate(all="ignore"):
        out = run(params, drivers, latitude_deg=float(config["site"]["latitude_deg"]),
                  frost_threshold_degc=frost, cold=cold, a7=a7)
    return out, sigma_m


def predictive_figure(out, data, config, label, n_lines: int = 50) -> None:
    observed = data[config["columns"]["nee"]].to_numpy(dtype=float)
    use = data[MASK_USE].to_numpy(dtype=bool)
    years = data.index.year.to_numpy()
    window = years < years[0] + 2
    days = data.index[window]
    fig, ax = plt.subplots(figsize=(FULL_WIDTH, 4.8))
    for line in out.nee[:n_lines, window]:
        ax.plot(days, line, color=BLUE, lw=0.6, alpha=0.18)
    ax.plot([], [], color=BLUE, lw=1, alpha=0.4, label=f"{n_lines} posterior runs of the model")
    ax.plot(days, np.median(out.nee[:, window], axis=0), color=BLUE, lw=1.8,
            label="posterior median")
    shown = window & use
    ax.scatter(data.index[shown], observed[shown], s=9, color=INK, zorder=3,
               label="observed NEE (days used in the fit)")
    ax.axhline(0, color=INK_SOFT, lw=0.8)
    ax.set_ylabel(f"NEE ({UNIT_FLUX})\nnegative = uptake")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncols=3)
    fig.tight_layout()
    save_figure(fig, f"fig_V_{label}_predictive_nee", config)


def pools_figure(out, data, config, label) -> None:
    first_year = int(data.index.year[0])
    fig, axes = plt.subplots(2, 3, figsize=(HALF_WIDTH, 5.4), sharex=True)
    for k, (ax, name) in enumerate(zip(axes.flat, POOLS)):
        values = out.pools_yearly[name]                      # (draws, years + 1)
        x = first_year + np.arange(values.shape[-1])
        low, mid, high = np.percentile(values, [5, 50, 95], axis=0)
        ax.fill_between(x, low, high, color=BLUE, alpha=0.2, lw=0)
        ax.plot(x, mid, color=BLUE, lw=2)
        ax.set_title(POOL_LABELS[name])
        ax.xaxis.set_major_locator(MaxNLocator(3, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        if k % 3 == 0:
            ax.set_ylabel("g C m$^{-2}$")
    handles = [plt.Line2D([], [], color=BLUE, lw=2, label="median"),
               plt.Line2D([], [], color=BLUE, lw=7, alpha=0.2, label="90% band")]
    fig.legend(handles=handles, loc="lower center", ncols=2, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save_figure(fig, f"fig_W_{label}_pools", config)


def uncertainty_figure(out, sigma_m, data, config, label) -> None:
    """RQ2: month by month, how large is each source of uncertainty? Parameter
    uncertainty is the spread of the model's NEE across posterior draws; measurement
    error is RANDUNC; model error is sigma_m (zero in the RANDUNC-only model)."""
    use = data[MASK_USE].to_numpy(dtype=bool)
    sigma = data[config["columns"]["sigma"]].to_numpy(dtype=float)
    months = data.index.month.to_numpy()
    parameter_sd = out.nee.std(axis=0)
    rows = []
    for m in range(1, 13):
        days = months == m
        rows.append((parameter_sd[days].mean(), np.nanmean(sigma[days & use]), float(np.median(sigma_m))))
    rows = np.array(rows)
    fig, ax = plt.subplots(figsize=(HALF_WIDTH, 5.0))
    x = np.arange(1, 13)
    series = (("parameter uncertainty", "parameters", BLUE, rows[:, 0]),
              ("measurement error (RANDUNC)", "RANDUNC", ORANGE, rows[:, 1]),
              ("model error ($\\sigma_m$)", "$\\sigma_m$", INK_SOFT, rows[:, 2]))
    for name, short, colour, values in series:
        ax.plot(x, values, color=colour, lw=2.4, marker="o", ms=6, label=name)
        ax.annotate(short, (12, values[-1]), xytext=(8, 0), textcoords="offset points",
                    va="center", fontsize=10.5, color=INK)
    ax.set_xticks(x, [m[0] for m in MONTHS])
    ax.set_xlim(0.5, 14.6)
    ax.set_ylim(bottom=0)
    ax.set_ylabel(f"standard deviation\n({UNIT_FLUX})")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncols=1)
    fig.tight_layout()
    save_figure(fig, f"fig_X_{label}_uncertainty", config)


# ---------------------------------------------------------------------------
# Summary table for slides
# ---------------------------------------------------------------------------

def latex_number(v: float) -> str:
    """Three significant figures: 0.0245, 137, 6,230; powers of ten below 0.001."""
    if v == 0 or not np.isfinite(v):
        return f"{v:g}"
    if abs(v) < 1e-3:
        mantissa, exponent = f"{v:.2e}".split("e")
        return rf"{float(mantissa):g}$\times$10\textsuperscript{{$-${abs(int(exponent))}}}"
    if abs(v) >= 1000:
        return f"{round(v, 2 - int(np.floor(np.log10(abs(v))))):,.0f}"
    return f"{v:.3g}"


def summary_table(idata, post, names, label, config) -> None:
    diagnostics = az.summary(idata, var_names=names)
    rows = []
    for name in names:
        values = chains_of(post, name).reshape(-1)
        rows.append({"quantity": name, "mean": values.mean(), "sd": values.std(),
                     "5%": np.percentile(values, 5), "95%": np.percentile(values, 95),
                     "r_hat": float(diagnostics.loc[name, "r_hat"]),
                     "ess_bulk": float(diagnostics.loc[name, "ess_bulk"]),
                     "ess_tail": float(diagnostics.loc[name, "ess_tail"])})
    table = pd.DataFrame(rows).set_index("quantity")
    out_dir = resolve_path(config["paths"]["figures"])
    out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_dir / f"{label}_summary.csv", float_format="%.6g")
    print(f"  saved {out_dir / f'{label}_summary.csv'}")

    # Slide version: fewer columns (the csv keeps them all). ESS is the smaller of
    # bulk and tail, as in the D49 rule; very small rates are shown in 10^k units.
    max_rhat = float(config["diagnostics"]["max_rhat"])
    min_ess = float(config["diagnostics"]["min_ess"])
    half = int(np.ceil(len(table) / 2))
    for part, chunk in enumerate((table.iloc[:half], table.iloc[half:]), start=1):
        lines = [r"\begin{tabular}[t]{lrrrrr}", r"\toprule",
                 r"parameter & mean & 5\% & 95\% & $\hat R$ & ESS \\", r"\midrule"]
        for name, r in chunk.iterrows():
            exponent = scale_exponent(chains_of(post, name)) if on_log_axis(name) else None
            scale = 1.0 if exponent is None else 10.0 ** exponent
            ess = min(r["ess_bulk"], r["ess_tail"])
            r_hat = f"{r['r_hat']:.3f}"
            lines.append(" & ".join([
                # \text{} (amsmath) sets subscripts in the slide's own font, as in the rest of the deck
                with_scale(short_label(name), exponent).replace("°", r"\textdegree{}")
                .replace(r"\mathrm{", r"\text{"),
                latex_number(r["mean"] / scale), latex_number(r["5%"] / scale),
                latex_number(r["95%"] / scale),
                r"\textbf{" + r_hat + "}" if r["r_hat"] >= max_rhat else r_hat,
                r"\textbf{" + f"{ess:,.0f}" + "}" if ess < min_ess else f"{ess:,.0f}",
            ]) + r" \\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        path = out_dir / f"{label}_summary_part{part}.tex"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"  saved {path}")


# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("label", help="name of the saved fit in results/, without .nc")
    parser.add_argument("--truth", default=None,
                        help="fake-data recovery table (csv with a 'truth' column) to mark true values")
    parser.add_argument("--predict", action="store_true",
                        help="real-data fits: also run the model forward (figures V, W, X)")
    parser.add_argument("--draws", type=int, default=200, help="posterior draws run forward")
    parser.add_argument("--per-page", type=int, default=6, help="parameters per trace page")
    args = parser.parse_args()

    config = load_config()
    apply_style()
    plt.rcParams.update(TALK_STYLE)
    idata = load(args.label, config)
    post = idata.posterior
    names = parameter_names(post)
    truth = read_truth(args.truth)
    heading(f"Talk figures for {args.label}: {len(names)} parameters, "
            f"{post[names[0]].shape[0]} chains x {post[names[0]].shape[1]} draws")

    summary_table(idata, post, names, args.label, config)
    trace_figures(post, names, args.label, config, truth, args.per_page)
    rank_figures(post, names, args.label, config)
    prior_posterior_figure(post, names, prior_draws(idata, names, config), args.label, config, truth)
    correlation_figures(post, names, args.label, config)

    if args.predict:
        if not args.label.startswith("fit"):
            print("  --predict skipped: only real-data fits (labels starting with 'fit') "
                  "have observed NEE to compare with")
            return
        data = fit_block(args.label, config)
        out, sigma_m = forward(idata, data, config, args.label, args.draws, seed=20261018)
        predictive_figure(out, data, config, args.label)
        pools_figure(out, data, config, args.label)
        uncertainty_figure(out, sigma_m, data, config, args.label)


if __name__ == "__main__":
    main()
