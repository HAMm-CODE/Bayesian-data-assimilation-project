"""Reduced models (D63): named sets of quantities held fixed at their prior medians.
The allocation is f_auto and three step fractions (D67), each fixed on its own.
Run with:  pytest"""

import numpy as np
import pytest
from scipy import stats

from dalec.config import load_config
from dalec.inference import SHARES, build_model
from dalec.priors import (
    FRACTIONS, apply_fixed, fixed_values, prior_median, priors_for, sample_prior, shares_from_fractions,
)
from test_inference import DRIVERS, LAT, fake_observations, hand_logp, model_logp, prior_logp
from test_model_numpy import median_parameters

M0 = ["c_woo_0", "c_som_0", "theta_som", "theta_woo", "frac_fol", "frac_lab", "frac_roo",
      "theta_min", "cr_onset", "lma", "f_auto"]


def test_config_fix_sets_are_valid_and_nested():
    sets = load_config()["fix_sets"]
    priors = priors_for(cold=True)
    for names in sets.values():
        fixed_values(names, priors)                       # raises on an unknown name
    assert set(sets["MA"]) < set(sets["MB"]) < set(sets["MC"]) < set(sets["M0"])
    assert sets["M0"] == M0


def test_fixed_values_are_prior_medians():
    priors = priors_for(cold=True)
    values = fixed_values(["c_woo_0", "f_auto", "frac_fol"], priors)
    assert values["c_woo_0"] == prior_median("c_woo_0", priors)
    assert values["f_auto"] == pytest.approx(stats.beta(25.5, 24.5).median())
    assert values["frac_fol"] == pytest.approx(stats.beta(13.06, 72.59).median())
    with pytest.raises(ValueError, match="cannot be fixed"):
        fixed_values(["f_lab"], priors)                    # a share of GPP is built, not sampled


def test_sample_prior_with_fixed_keeps_the_other_draws():
    """Fixing is applied after the random draws, so the free quantities are unchanged."""
    priors = priors_for(cold=True)
    fixed = fixed_values(M0, priors)
    _, raw_free = sample_prior(50, np.random.default_rng(3), priors)
    params, raw = sample_prior(50, np.random.default_rng(3), priors, fixed=fixed)
    for name in ("ceff", "c_lf", "d_onset", "c_fol_0", "theta_lit"):
        assert np.array_equal(raw[name], raw_free[name])
    for name in M0:
        assert np.all(raw[name] == fixed[name])
    expected = shares_from_fractions(*(fixed[k] for k in ("f_auto", *FRACTIONS)))
    for name in SHARES:
        assert np.allclose(raw[name], expected[name])
        assert np.allclose(getattr(params, name), raw[name])
    assert np.allclose(params.theta_min, fixed["theta_min"])


def test_fixing_f_auto_alone_keeps_each_draws_fractions():
    priors = priors_for(cold=False)
    _, raw = sample_prior(20, np.random.default_rng(4), priors)
    _, fixed_raw = apply_fixed(raw, {"f_auto": 0.5})
    for name in FRACTIONS:
        assert np.allclose(raw[name], fixed_raw[name])
    total = sum(fixed_raw[n] for n in SHARES)
    assert np.allclose(total, 0.5)


def test_m0_log_posterior_equals_hand_computation_minus_fixed_priors():
    """D63: the reduced model's log density is the full one at the same point, minus the
    priors of the fixed quantities."""
    priors = priors_for(cold=True)
    fixed = fixed_values(M0, priors)
    shares = shares_from_fractions(*(fixed[k] for k in ("f_auto", *FRACTIONS)))
    params = median_parameters(
        f_auto=fixed["f_auto"], lma=fixed["lma"], cr_onset=fixed["cr_onset"],
        theta_som=fixed["theta_som"], theta_woo=fixed["theta_woo"], theta_min=fixed["theta_min"],
        c_woo_0=fixed["c_woo_0"], c_som_0=fixed["c_som_0"], **{k: float(v) for k, v in shares.items()})
    observed, sigma, use = fake_observations(params)
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True, fixed=fixed)
    free_names = {rv.name for rv in model.free_RVs}
    assert not free_names & set(M0)
    assert len(free_names) == 14 + 1                       # 14 quantities + sigma_m (D63)

    full, point = hand_logp(params, 0.4, observed, sigma, use,
                            "randunc_plus_model_error", cold=(-6.0, 5.0))
    removed = sum(prior_logp(name, fixed[name], priors) for name in M0)
    free_point = {k: v for k, v in point.items() if k in free_names}
    assert model_logp(model, free_point) == pytest.approx(full - removed, rel=1e-9)


def test_allowed_misses_matches_d49_d59_and_gives_3_for_m0():
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / "scripts" / "04_fake_data_check.py"
    spec = importlib.util.spec_from_file_location("fake_check_script", path)
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    tail = load_config()["fake_data"]["miss_tail"]
    got = {k: script.allowed_misses(k, 0.1, tail) for k in (14, 22, 24, 25, 26)}
    assert got == {14: 3, 22: 4, 24: 4, 25: 4, 26: 4}
    assert 24 - got[24] == load_config()["fake_data"]["min_covered"]     # same as D49
