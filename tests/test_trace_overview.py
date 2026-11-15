"""Trace overview (script 15): one panel per quantity, thesis labels, no code names. Run with:  pytest"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import xarray as xr

_path = Path(__file__).resolve().parents[1] / "scripts" / "15_trace_overview.py"
_spec = importlib.util.spec_from_file_location("trace_overview_script", _path)
script = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(script)


def posterior(names, chains=4, draws=50):
    rng = np.random.default_rng(1)
    return xr.Dataset({n: (("chain", "draw"), 1e-4 + rng.random((chains, draws)) * 1e-4 if n.startswith("theta")
                           else rng.random((chains, draws)) + 1) for n in names})


def test_one_panel_per_quantity_and_spare_panels_hidden():
    names = ["f_auto", "frac_fol", "theta_roo", "lma", "sigma_m", "cr_onset", "c_lab_0"]
    fig = script.trace_overview_figure(posterior(names), names, columns=4)
    used = [ax for ax in fig.axes if ax.axison]
    assert len(used) == len(names)
    assert len(fig.axes) == 8                                # 2 rows x 4 columns, one spare panel


def test_titles_use_the_notation_sheet_not_code_names():
    names = ["frac_fol", "lma", "cr_onset", "c_lab_0"]
    fig = script.trace_overview_figure(posterior(names), names, columns=2)
    titles = " ".join(ax.get_title(loc="left") for ax in fig.axes if ax.axison)
    assert "s_" not in titles and "LMA" not in titles and "cr_" not in titles
    assert r"c_\mathrm{ronset}" in titles and r"c_\mathrm{lma}" in titles and r"C^0_\mathrm{lab}" in titles


def test_all_22_fit_on_a_6_by_4_grid():
    names = [f"c_{i}" for i in range(22)]
    fig = script.trace_overview_figure(posterior(names), names)
    assert len(fig.axes) == 24
