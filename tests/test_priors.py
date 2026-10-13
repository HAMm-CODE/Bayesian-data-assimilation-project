"""Tests for src/dalec/priors.py. Run from the repository root with:  pytest"""

import numpy as np
import pytest

from dalec.priors import C_LF_RANGE, PRIORS, sample_prior

N = 20000


@pytest.fixture(scope="module")
def draws():
    return sample_prior(N, np.random.default_rng(0))


def test_allocation_adds_up_to_one(draws):
    params, _ = draws
    total = params.f_auto + params.f_lab + params.f_fol + params.f_roo + params.f_woo
    assert np.allclose(total, 1.0)


def test_lognormal_medians_match_decisions(draws):
    """Each log-normal prior's sample median is close to the median written in D35 to D39."""
    _, raw = draws
    for name, (kind, median, _) in PRIORS.items():
        if kind == "lognormal":
            assert np.median(raw[name]) == pytest.approx(median, rel=0.03), name


def test_beta_and_normal_means_match_decisions(draws):
    _, raw = draws
    assert raw["f_auto"].mean() == pytest.approx(0.51, abs=0.005)
    assert raw["h"].mean() == pytest.approx(0.20, abs=0.005)
    assert raw["d_onset"].mean() == pytest.approx(145, abs=0.5)
    assert raw["d_fall"].mean() == pytest.approx(265, abs=0.5)


def test_theta_min_follows_from_h(draws):
    _, raw = draws
    h_back = raw["theta_min"] / (raw["theta_lit"] + raw["theta_min"])
    assert np.allclose(h_back, raw["h"])


def test_c_lf_stays_inside_polynomial_range(draws):
    params, _ = draws
    low, high = C_LF_RANGE
    assert np.all(params.c_lf >= low) and np.all(params.c_lf < high)


def test_same_seed_gives_same_draws():
    a, _ = sample_prior(50, np.random.default_rng(7))
    b, _ = sample_prior(50, np.random.default_rng(7))
    assert np.array_equal(a.ceff, b.ceff) and np.array_equal(a.c_som_0, b.c_som_0)


# --- soft steady-state constraint on the foliage (D44) -----------------------

from dalec.priors import EDC_EQF, edc_log_weight, steady_state_ratio  # noqa: E402


def test_balanced_pool_has_no_penalty():
    assert edc_log_weight(steady_state_ratio(100.0, 100.0, 50.0, 50.0)) == 0.0


def test_penalty_is_symmetric_and_gaussian_in_log_rs():
    """A factor EQF too high or too low costs -0.5; two factors cost -2 (EDC2_400.c EDCPROB)."""
    assert edc_log_weight(EDC_EQF) == pytest.approx(-0.5)
    assert edc_log_weight(1 / EDC_EQF) == pytest.approx(-0.5)
    assert edc_log_weight(EDC_EQF**2) == pytest.approx(-2.0)


def test_growing_and_shrinking_pools_are_penalised():
    growing = steady_state_ratio(300.0, 100.0, 80.0, 50.0)   # inputs exceed outputs, pool rising
    shrinking = steady_state_ratio(50.0, 200.0, 20.0, 50.0)
    assert growing > 1 and edc_log_weight(growing) < -0.5
    assert shrinking < 1 and edc_log_weight(shrinking) < -0.5


def test_broken_runs_get_zero_weight():
    assert edc_log_weight(steady_state_ratio(100.0, 0.0, 50.0, 50.0) * -1) == -np.inf
    assert edc_log_weight(np.nan) == -np.inf