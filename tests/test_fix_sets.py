"""Reduced models (D63): named sets of quantities held fixed at their prior medians,
with the four growth shares fixed together. Run with:  pytest"""

import numpy as np
import pytest
from scipy import stats

from dalec.config import load_config
from dalec.inference import SHARES, build_model
from dalec.priors import (
    DIRICHLET_ALPHA, SHARES_KEY, apply_fixed, fixed_values, prior_median, priors_for, sample_prior,
    shares_mean,
)
from test_inference import DRIVERS, LAT, fake_observations, hand_logp, model_logp, prior_logp
from test_model_numpy import median_parameters

M0 = ["c_woo_0", "c_som_0", "theta_som", "theta_woo", "shares", "h", "cr_onset", "lma", "f_auto"]


def test_config_fix_sets_are_valid_and_nested():
    sets = load_config()["fix_sets"]
    priors = priors_for(cold=True)
    for names in sets.values():
        fixed_values(names, priors)                       # raises on an unknown name
    assert set(sets["MA"]) < set(sets["MB"]) < set(sets["MC"]) < set(sets["M0"])
    assert sets["M0"] == M0


def test_fixed_values_are_prior_medians_and_mean_shares():
    priors = priors_for(cold=True)
    values = fixed_values(["c_woo_0", "f_auto", SHARES_KEY], priors)
    assert values["c_woo_0"] == prior_median("c_woo_0", priors)
    assert values["f_auto"] == pytest.approx(stats.beta(25.5, 24.5).median())
    assert np.allclose(values[SHARES_KEY], shares_mean()) and shares_mean().sum() == pytest.approx(1.0)
    with pytest.raises(ValueError, match="cannot be fixed"):
        fixed_values(["f_lab"], priors)                    # a single share cannot be fixed alone


def test_sample_prior_with_fixed_keeps_the_other_draws():
    """Fixing is applied after the random draws, so the free quantities are unchanged."""
    priors = priors_for(cold=True)
    fixed = fixed_values(M0, priors)
    _, raw_free = sample_prior(50, np.random.default_rng(3), priors)
    params, raw = sample_prior(50, np.random.default_rng(3), priors, fixed=fixed)
    for name in ("ceff", "c_lf", "d_onset", "c_fol_0", "theta_lit"):
        assert np.array_equal(raw[name], raw_free[name])
    for name in ("c_woo_0", "lma", "f_auto", "h"):
        assert np.all(raw[name] == fixed[name])
    growth = 1.0 - raw["f_auto"]
    for k, name in enumerate(SHARES):
        assert np.allclose(raw[name], growth * shares_mean()[k])
        assert np.allclose(getattr(params, name), raw[name])
    assert np.allclose(raw["theta_min"], raw["theta_lit"] * raw["h"] / (1 - raw["h"]))
    assert np.allclose(params.theta_min, raw["theta_min"])


def test_fixing_f_auto_alone_keeps_each_draws_shares():
    priors = priors_for(cold=False)
    _, raw = sample_prior(20, np.random.default_rng(4), priors)
    _, fixed_raw = apply_fixed(raw, {"f_auto": 0.5})
    old = np.stack([raw[n] for n in DIRICHLET_ALPHA], 1) / (1 - raw["f_auto"])[:, None]
    new = np.stack([fixed_raw[n] for n in DIRICHLET_ALPHA], 1) / 0.5
    assert np.allclose(old, new)


def test_m0_log_posterior_equals_hand_computation_minus_fixed_priors():
    """D63: the reduced model's log density is the full one at the same point, minus the
    priors of the fixed quantities and minus the Dirichlet term of the fixed shares."""
    priors = priors_for(cold=True)
    fixed = fixed_values(M0, priors)
    growth = 1.0 - fixed["f_auto"]
    params = median_parameters(
        f_auto=fixed["f_auto"], lma=fixed["lma"], cr_onset=fixed["cr_onset"],
        theta_som=fixed["theta_som"], theta_woo=fixed["theta_woo"],
        c_woo_0=fixed["c_woo_0"], c_som_0=fixed["c_som_0"],
        **{name: growth * fixed[SHARES_KEY][k] for k, name in enumerate(SHARES)})
    observed, sigma, use = fake_observations(params)
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True, fixed=fixed)
    free_names = {rv.name for rv in model.free_RVs}
    assert not free_names & set(M0)
    assert SHARES_KEY not in free_names
    assert len(free_names) == 14 + 1                       # 14 quantities + sigma_m (D63)

    full, point = hand_logp(params, fixed["h"], 0.4, observed, sigma, use,
                            "randunc_plus_model_error", cold=(-6.0, 5.0))
    removed = sum(prior_logp(name, fixed[name], priors) for name in M0 if name != SHARES_KEY)
    removed += stats.dirichlet(list(DIRICHLET_ALPHA.values())).logpdf(fixed[SHARES_KEY])
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
