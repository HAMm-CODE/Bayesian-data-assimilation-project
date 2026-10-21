"""Checks of the scoring helpers in scripts/06_evaluate.py (D61)."""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

_path = Path(__file__).resolve().parents[1] / "scripts" / "06_evaluate.py"
_spec = importlib.util.spec_from_file_location("evaluate_script", _path)
ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ev)


def test_lpd_is_log_of_average_density():
    rng = np.random.default_rng(0)
    loglik = rng.normal(-1.0, 0.5, size=(50, 7))
    expected = np.log(np.exp(loglik).mean(axis=0))
    assert np.allclose(ev.lpd_per_day(loglik), expected)


def test_pointwise_loglik_uses_randunc_and_model_error():
    model = np.array([[0.0, 1.0], [0.5, 1.5]])          # 2 draws x 2 days
    observed = np.array([0.2, 1.1])
    sigma = np.array([0.3, 0.4])
    sigma_m = np.array([0.0, 0.5])
    got = ev.pointwise_loglik(model, observed, sigma, sigma_m)
    sd = np.sqrt(0.4 ** 2 + 0.5 ** 2)
    hand = -0.5 * np.log(2 * np.pi * sd ** 2) - 0.5 * ((1.1 - 1.5) / sd) ** 2
    assert np.isclose(got[1, 1], hand)
    assert np.isclose(got[0, 0], -0.5 * np.log(2 * np.pi * 0.09) - 0.5 * (0.2 / 0.3) ** 2)


def test_widened_se_grows_with_positive_correlation_and_skips_gaps():
    rng = np.random.default_rng(1)
    n, rho = 4000, 0.5
    x = np.empty(n)
    x[0] = rng.normal()
    for t in range(1, n):
        x[t] = rho * x[t - 1] + rng.normal() * np.sqrt(1 - rho ** 2)
    dates = pd.date_range("2011-01-01", periods=n, freq="D")
    plain, r, wide = ev.widened_se(x, dates)
    assert abs(r - rho) < 0.05
    assert np.isclose(wide / plain, np.sqrt((1 + r) / (1 - r)))
    # with every other day missing there are no consecutive pairs left: no widening
    plain2, r2, wide2 = ev.widened_se(x[::2], dates[::2])
    assert r2 == 0.0 and wide2 == plain2


def test_negative_correlation_does_not_shrink_the_error():
    diff = np.array([1.0, -1.0] * 50)
    dates = pd.date_range("2011-01-01", periods=diff.size, freq="D")
    plain, r, wide = ev.widened_se(diff, dates)
    assert r < 0 and wide == plain
