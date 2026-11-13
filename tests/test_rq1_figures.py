"""RQ1 figures (script 13): the numbers behind the two figures are the ones in
the saved tables, and both figures can be drawn. Run with:  pytest"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd
import pytest

_path = Path(__file__).resolve().parents[1] / "scripts" / "13_rq1_figures.py"
_spec = importlib.util.spec_from_file_location("rq1_figures_script", _path)
script = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(script)


def parameters() -> pd.DataFrame:
    return pd.DataFrame(
        {"width_ratio": [0.16, 0.50, 1.04, 0.0, 0.54],
         "narrowing": ["strongly informed", "partly informed", "not informed", "fixed", "partly informed"]},
        index=pd.Index(["ceff", "f_auto", "c_roo_0", "c_woo_0", "f_woo"], name="quantity"))


def own() -> pd.DataFrame:
    return pd.DataFrame(
        {"own_prior": [0.19, 0.01, 0.13], "likelihood": [0.28, 0.31, 0.02],
         "reading": ["own prior and data pull against each other", "data decide; own prior hardly matters",
                     "own prior decides, data weak"]},
        index=["ceff", "cold_t_full_fraction", "f_roo"])


def test_reduction_is_one_minus_width_ratio_in_percent():
    reduction = script.reduction_table(parameters())
    assert reduction["ceff"] == pytest.approx(84.0)
    assert reduction["f_auto"] == pytest.approx(50.0)
    assert reduction["c_roo_0"] == pytest.approx(-4.0)          # a posterior wider than its prior gives a negative bar


def test_fixed_and_derived_quantities_get_no_bar():
    reduction = script.reduction_table(parameters())
    assert "c_woo_0" not in reduction.index                      # fixed in the model
    assert "f_woo" not in reduction.index                        # computed, not sampled
    assert list(reduction.index) == ["c_roo_0", "f_auto", "ceff"]   # smallest first, so the largest is drawn on top


def test_sensitivity_table_keeps_the_saved_numbers_and_renames_cold_thresholds():
    table = script.sensitivity_table(own())
    assert list(table.index) == ["cold_t_full", "f_roo", "ceff"]    # sorted by own-prior sensitivity
    assert table.loc["ceff", "own_prior"] == pytest.approx(0.19)
    assert table.loc["f_roo", "likelihood"] == pytest.approx(0.02)


def test_both_figures_draw():
    assert len(script.reduction_figure(parameters()).axes[0].patches) == 3
    assert len(script.sensitivity_figure(own()).axes[0].collections) == 2
