"""Prior sets (D64): a named set of changed priors, used everywhere the site
priors are used: draws, fixed values, the PyMC model and the saved fits.
These tests use a made-up set so they do not depend on the literature values.
Run with:  pytest"""

import numpy as np
import pytest
from scipy import stats

import dalec.priors as priors_module
from dalec.inference import build_model
from dalec.priors import (
    FRACTIONS, fixed_values, log_prior, priors_for, priors_from_attrs, priors_to_text, sample_prior,
)
from test_inference import DRIVERS, LAT, fake_observations, hand_logp, model_logp, prior_logp
from test_model_numpy import median_parameters

TEST_SET = {
    "c_woo_0": ("lognormal", 8000.0, 0.5),
    "theta_som": ("lognormal", 4.0e-5, 0.6),
    "f_auto": ("beta", 12.0, 12.0),
    "lma": ("lognormal", 120.0, 0.4),
    "frac_fol": ("beta", 5.0, 30.0),
}


@pytest.fixture
def test_set(monkeypatch):
    monkeypatch.setitem(priors_module.PRIOR_SETS, "test", dict(TEST_SET))
    return "test"


def test_unknown_sets_are_refused():
    with pytest.raises(ValueError, match="unknown prior set"):
        priors_for(True, "nonsense")
    assert priors_for(True, "site") == priors_for(True)


def test_set_changes_only_its_quantities(test_set):
    site, changed = priors_for(True), priors_for(True, test_set)
    assert set(changed) == set(site)
    for name in site:
        assert changed[name] == TEST_SET.get(name, site[name])


def test_draws_and_fixed_values_follow_the_set(test_set):
    priors = priors_for(True, test_set)
    _, raw = sample_prior(40_000, np.random.default_rng(5), priors)
    assert np.median(raw["c_woo_0"]) == pytest.approx(8000.0, rel=0.02)
    assert raw["frac_fol"].mean() == pytest.approx(5.0 / 35.0, abs=0.003)
    values = fixed_values(["c_woo_0", "f_auto", "frac_fol"], priors)
    assert values["c_woo_0"] == 8000.0
    assert values["f_auto"] == pytest.approx(0.5)
    assert values["frac_fol"] == pytest.approx(stats.beta(5.0, 30.0).median())


def test_set_survives_saving_with_a_fit(test_set):
    priors = priors_for(True, test_set)
    back = priors_from_attrs({"priors": priors_to_text(priors)})
    assert back == priors
    _, raw = sample_prior(10, np.random.default_rng(6), back)
    assert np.isfinite(log_prior(raw, back)).all()


def test_log_prior_uses_the_set(test_set):
    priors = priors_for(True, test_set)
    _, raw = sample_prior(3, np.random.default_rng(7), priors)
    by_hand = np.array([sum(prior_logp(name, raw[name][i], priors) for name in priors) for i in range(3)])
    assert np.allclose(log_prior(raw, priors), by_hand, rtol=1e-10)


def test_model_log_density_changes_by_exactly_the_prior_difference(test_set):
    """Same point, same data: the model with the test set differs from the site model
    only by the prior terms of the changed quantities."""
    params = median_parameters()
    observed, sigma, use = fake_observations(params)
    site = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                       observation_model="randunc_plus_model_error", cold_limit=True)
    lit = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                      observation_model="randunc_plus_model_error", cold_limit=True,
                      prior_set=test_set)
    _, point = hand_logp(params, 0.4, observed, sigma, use, "randunc_plus_model_error",
                         cold=(-6.0, 5.0))
    p_site, p_lit = priors_for(True), priors_for(True, test_set)
    expected = sum(prior_logp(n, point[n], p_lit) - prior_logp(n, point[n], p_site) for n in TEST_SET)
    assert model_logp(lit, point) - model_logp(site, point) == pytest.approx(expected, rel=1e-9)


def test_reduced_model_with_a_set(test_set):
    """A fix set with a prior set: fixed at the set's medians (D64, L models)."""
    priors = priors_for(True, test_set)
    fixed = fixed_values(["c_woo_0", *FRACTIONS], priors)
    observed, sigma, use = fake_observations(median_parameters())
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True,
                        fixed=fixed, prior_set=test_set)
    free = {rv.name for rv in model.free_RVs}
    assert "c_woo_0" not in free and not free & set(FRACTIONS) and "lma" in free
    assert np.isfinite(model.compile_logp()(model.initial_point()))


def test_the_literature_sets_are_valid():
    """The filled-in sets (D64, F35, D67) load, change only known quantities, and the
    exploratory sets change only the two starting pools."""
    for name in ("lit", "lit_wood", "lit_wood_soil"):
        assert set(priors_for(True, name)) == set(priors_for(True))
    assert set(priors_module.PRIOR_SETS["lit_wood_soil"]) == {"c_woo_0", "c_som_0"}
    lit = priors_module.PRIOR_SETS["lit"]
    assert lit["c_woo_0"] == priors_module.PRIOR_SETS["lit_wood"]["c_woo_0"]
    values = fixed_values(["c_woo_0", "c_som_0"], priors_for(True, "lit_wood_soil"))
    assert values == {"c_woo_0": lit["c_woo_0"][1], "c_som_0": lit["c_som_0"][1]}
