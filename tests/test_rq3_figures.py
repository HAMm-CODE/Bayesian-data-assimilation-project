"""RQ3 figures (script 14): the recovery plot uses the saved table on the prior's
scale, the pairwise plot picks the strongly tied quantities, and both draw.
Run with:  pytest"""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import pytest

_path = Path(__file__).resolve().parents[1] / "scripts" / "14_rq3_figures.py"
_spec = importlib.util.spec_from_file_location("rq3_figures_script", _path)
script = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(script)


def recovery() -> pd.DataFrame:
    rows = {
        # a log-normal quantity whose posterior equals its prior and whose truth is the prior's 95% point
        "ceff": dict(truth=40.0, prior_low=10.0, prior_median=20.0, prior_high=40.0, posterior_low=10.0,
                     posterior_median=20.0, posterior_high=40.0, truth_inside=True, width_ratio=1.0,
                     narrowing="not informed"),
        "c_woo_0": dict(truth=5000.0, prior_low=5000.0, prior_median=5000.0, prior_high=5000.0, posterior_low=5000.0,
                        posterior_median=5000.0, posterior_high=5000.0, truth_inside=True, width_ratio=0.0,
                        narrowing="fixed"),
        "f_auto": dict(truth=0.7, prior_low=0.4, prior_median=0.5, prior_high=0.6, posterior_low=0.45,
                       posterior_median=0.5, posterior_high=0.55, truth_inside=False, width_ratio=0.5,
                       narrowing="partly informed"),
    }
    return pd.DataFrame(rows).T.rename_axis("quantity")


def test_prior_units_put_the_prior_range_at_minus_one_to_one():
    units = script.in_prior_units(recovery())
    assert units.loc["ceff", "median"] == pytest.approx(0.0)
    assert units.loc["ceff", "low"] == pytest.approx(-1.0)        # log scale: 10, 20, 40 are evenly spaced
    assert units.loc["ceff", "truth"] == pytest.approx(1.0)
    assert units.loc["f_auto", "truth"] > 1.0                      # outside the prior's 90% range
    assert "c_woo_0" not in units.index                            # fixed quantities get no row


def test_tied_quantities_are_the_ones_with_a_strong_partner():
    names = ["a", "b", "c", "d"]
    corr = pd.DataFrame(np.eye(4), index=names, columns=names)
    corr.loc["a", "b"] = corr.loc["b", "a"] = -0.8
    corr.loc["c", "d"] = corr.loc["d", "c"] = 0.3
    assert script.tied_quantities(corr) == ["a", "b"]
    assert script.tied_quantities(corr, most=1)                    # never fewer than two: falls back to the most tied


def test_both_figures_draw():
    table = recovery()
    fig = script.recovery_figure([table, table], ["true forest 1", "true forest 2"])
    assert len(fig.axes) == 2
    rng = np.random.default_rng(0)
    z = rng.multivariate_normal([0, 0, 0], [[1, -0.8, 0], [-0.8, 1, 0], [0, 0, 1]], 500)
    draws = {"ceff": np.exp(z[:, 0]), "lma": np.exp(z[:, 1]), "f_auto": 1 / (1 + np.exp(-z[:, 2]))}
    corr = pd.DataFrame(z, columns=list(draws)).corr()
    fig = script.pairwise_figure(draws, corr)
    assert len(fig.axes) == 9
    # correlations only: no true-value marker is ever drawn in the pairwise plot
    from matplotlib.colors import to_hex
    orange = to_hex(script.ORANGE)
    for ax in fig.axes:
        assert all(to_hex(line.get_color()) != orange for line in ax.lines)
        for cloud in ax.collections:
            assert all(to_hex(c) != orange for c in cloud.get_facecolor())


def test_recovery_density_figure_marks_truth_and_misses():
    table = recovery()
    rng = np.random.default_rng(1)
    prior = {"ceff": np.exp(rng.normal(3.0, 0.4, 2000)), "f_auto": rng.beta(20, 20, 2000)}
    posterior = {"ceff": np.exp(rng.normal(3.0, 0.2, 2000)), "f_auto": rng.beta(80, 80, 2000)}
    fig = script.recovery_density_figure(prior, posterior, table)
    drawn = [ax for ax in fig.axes if ax.get_title()]
    assert len(drawn) == 2                                           # the fixed quantity has no panel
    texts = [" ".join(t.get_text() for t in ax.texts) for ax in drawn]
    assert sum("missed" in text for text in texts) == 1             # f_auto: truth outside the posterior range
