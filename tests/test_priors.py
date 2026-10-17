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


# ---------------------------------------------------------------------------
# The cold-limit model's priors (D55, D57)
# ---------------------------------------------------------------------------

from dalec.inference import initial_values  # noqa: E402
from dalec.priors import (  # noqa: E402
    COLD_MODEL_PRIOR_CHANGES, COLD_PRIORS, log_prior, priors_for, priors_from_attrs,
    priors_to_text, sample_cold, to_prior_scale,
)


def test_cold_model_changes_only_ceff_to_the_d38_prior():
    cold, main = priors_for(cold=True), priors_for(cold=False)
    assert main == PRIORS
    assert cold["ceff"] == ("lognormal", 16.5, 0.41)
    assert {k for k in PRIORS if cold[k] != main[k]} == set(COLD_MODEL_PRIOR_CHANGES) == {"ceff"}


def test_saved_prior_list_round_trips_and_old_fits_mean_main_priors():
    cold = priors_for(cold=True)
    assert priors_from_attrs({"priors": priors_to_text(cold)}) == cold
    assert priors_from_attrs({}) == PRIORS      # fits saved before D57 used PRIORS


def test_cold_thresholds_stay_in_their_ranges_and_centre_on_zero_scale():
    draws = sample_cold(5000, np.random.default_rng(0))
    for name, (low, high, _, _) in COLD_PRIORS.items():
        assert low < draws[name].min() and draws[name].max() < high
        assert abs(float(to_prior_scale(name, (low + high) / 2))) < 1e-9


def test_log_prior_uses_the_given_prior_list():
    _, raw = sample_prior(50, np.random.default_rng(1), priors_for(cold=True))
    from scipy import stats
    diff = log_prior(raw, priors_for(cold=True)) - log_prior(raw)
    expected = (stats.lognorm(s=0.41, scale=16.5).logpdf(raw["ceff"])
                - stats.lognorm(s=0.354, scale=10.0).logpdf(raw["ceff"]))
    assert np.allclose(diff, expected)


def test_chain_start_includes_cold_fractions():
    rng = np.random.default_rng(2)
    _, raw = sample_prior(3, rng, priors_for(cold=True))
    raw.update(sample_cold(3, rng))
    start = initial_values(raw, 1)
    for name, (low, high, _, _) in COLD_PRIORS.items():
        assert start[f"{name}_fraction"] == pytest.approx((raw[name][1] - low) / (high - low))