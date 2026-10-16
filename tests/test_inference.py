"""The PyMC model must add up exactly to: priors + soft constraint + likelihood,
computed independently with scipy and the numpy model. Run with:  pytest"""

import numpy as np
import pytest
from scipy import stats

from dalec.inference import SHARES, build_model, initial_values
from dalec.model_numpy import run
from dalec.priors import (
    C_LF_RANGE, COLD_PRIORS, DIRICHLET_ALPHA, PRIORS, edc_log_weight, log_prior, sample_prior,
    steady_state_ratio,
)
from test_model_numpy import LAT, median_parameters, synthetic_drivers

DRIVERS = synthetic_drivers(2)


def prior_logp(name, value):
    kind, a, b = PRIORS[name]
    if kind == "lognormal":
        return stats.lognorm(s=b, scale=a).logpdf(value)
    if kind == "normal":
        return stats.norm(a, b).logpdf(value)
    dist = stats.beta(a, b)
    if name == "c_lf":
        low, high = C_LF_RANGE
        return dist.logpdf(value) - np.log(dist.cdf(high) - dist.cdf(low))
    return dist.logpdf(value)


def fake_observations(params, seed=1):
    rng = np.random.default_rng(seed)
    nee = run(params, DRIVERS, latitude_deg=LAT).nee
    sigma = 0.3 + 0.2 * rng.random(nee.size)
    use = rng.random(nee.size) > 0.2
    observed = np.where(use, nee + sigma * rng.standard_normal(nee.size), np.nan)
    return observed, sigma, use


def hand_logp(params, h, sigma_m, observed, sigma, use, observation_model, cold=None):
    """Priors + D44 term + likelihood, computed with scipy and the numpy model.
    `cold` = (t_shutdown, gap) adds the D55 cold limit and its two priors."""
    values = {name: float(getattr(params, name)) for name in PRIORS if name != "h"}
    values["h"] = h
    growth = 1.0 - values["f_auto"]
    shares = np.array([float(getattr(params, s)) / growth for s in SHARES])
    # theta_min in the model follows from h; make the numpy run use the same value.
    theta_min = values["theta_lit"] * h / (1 - h)
    reference = run(type(params)(**{**params.as_arrays(), "theta_min": theta_min}),
                    DRIVERS, latitude_deg=LAT, cold=cold)

    total = sum(prior_logp(name, value) for name, value in values.items())
    if cold is not None:
        # The model samples each threshold as a Beta(2, 2) fraction of its range.
        for value, (low, high, a, b) in zip(cold, COLD_PRIORS.values()):
            total += stats.beta(a, b).logpdf((value - low) / (high - low))
    total += stats.dirichlet(list(DIRICHLET_ALPHA.values())).logpdf(shares)
    rs = steady_state_ratio(reference.pools_in["c_fol"], reference.pools_out["c_fol"],
                            reference.pools_yearly["c_fol"].mean(), reference.pools_start["c_fol"])
    total += -0.5 * (np.log(rs) / np.log(2.0)) ** 2
    sd = sigma[use]
    if observation_model == "randunc_plus_model_error":
        total += stats.halfnorm(scale=1.0).logpdf(sigma_m)
        sd = np.sqrt(sd**2 + sigma_m**2)
    total += stats.norm(reference.nee[use], sd).logpdf(observed[use]).sum()
    point = {**values, "shares": shares, "sigma_m": sigma_m}
    if cold is not None:
        for value, (name, (low, high, _, _)) in zip(cold, COLD_PRIORS.items()):
            point[f"{name}_fraction"] = (value - low) / (high - low)
    return float(total), point


def model_logp(model, point):
    """PyMC's log density at `point`, without the change-of-scale (Jacobian)
    terms, so it is the plain density on the parameters' own scale."""
    transformed = {}
    for rv in model.free_RVs:
        # float64 on purpose: PyTensor stores a plain Python number such as 16.5,
        # which fits exactly in 32 bits, as a 32-bit number, and then rounds the
        # result of the log to about 7 digits.
        value = np.asarray(point[rv.name], dtype=np.float64)
        transform = model.rvs_to_transforms[rv]
        transformed[model.rvs_to_values[rv].name] = (
            value if transform is None else transform.forward(value, *rv.owner.inputs).eval())
    return float(model.compile_logp(jacobian=False)(transformed))


@pytest.mark.parametrize("observation_model", ["randunc_only", "randunc_plus_model_error"])
def test_log_posterior_equals_hand_computation(observation_model):
    observed, sigma, use = fake_observations(median_parameters())
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model=observation_model)
    a_expected, a_point = hand_logp(median_parameters(), 0.2, 0.4, observed, sigma, use,
                                    observation_model)
    other = median_parameters(ceff=12.0, lma=160.0, c_lf=0.3, temperature_exponent=0.04,
                              d_fall=255.0, c_som_0=6500.0)
    b_expected, b_point = hand_logp(other, 0.25, 0.7, observed, sigma, use, observation_model)
    a, b = model_logp(model, a_point), model_logp(model, b_point)
    assert a == pytest.approx(a_expected, rel=1e-9)
    assert b == pytest.approx(b_expected, rel=1e-9)


def test_initial_values_round_trip():
    _, raw = sample_prior(5, np.random.default_rng(3))
    start = initial_values(raw, 2)
    assert start["shares"].sum() == pytest.approx(1.0)
    growth = 1 - start["f_auto"]
    for k, name in enumerate(SHARES):
        assert growth * start["shares"][k] == pytest.approx(raw[name][2])


def test_rejects_assimilable_day_without_sigma():
    observed, sigma, use = fake_observations(median_parameters())
    sigma[np.flatnonzero(use)[0]] = 0.0
    with pytest.raises(ValueError, match="sigma"):
        build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT)


def test_log_prior_plus_terms_equals_hand_computation():
    """priors.log_prior (used by the chain diagnosis) + D44 term + likelihood
    equals the independent hand computation of the full log posterior."""
    params = median_parameters()
    observed, sigma, use = fake_observations(params)
    expected, point = hand_logp(params, 0.2, 0.0, observed, sigma, use, "randunc_only")
    raw = {name: point[name] for name in PRIORS}
    growth = 1.0 - point["f_auto"]
    for k, name in enumerate(SHARES):
        raw[name] = growth * point["shares"][k]
    reference = run(type(params)(**{**params.as_arrays(),
                                    "theta_min": point["theta_lit"] * 0.2 / 0.8}),
                    DRIVERS, latitude_deg=LAT)
    rs = steady_state_ratio(reference.pools_in["c_fol"], reference.pools_out["c_fol"],
                            reference.pools_yearly["c_fol"].mean(), reference.pools_start["c_fol"])
    likelihood = stats.norm(reference.nee[use], sigma[use]).logpdf(observed[use]).sum()
    total = float(log_prior(raw)) + float(edc_log_weight(rs)) + likelihood
    assert total == pytest.approx(expected, rel=1e-9)


def test_fixed_parameters_are_constants_and_drop_out_of_the_prior():
    """D54: a fixed parameter is not sampled, and the log density equals the full
    one minus that parameter's own prior term."""
    params = median_parameters(lma=185.0)
    observed, sigma, use = fake_observations(params)
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT, fixed={"lma": 185.0})
    assert "lma" not in {rv.name for rv in model.free_RVs}
    expected, point = hand_logp(params, 0.2, 0.0, observed, sigma, use, "randunc_only")
    free = {k: v for k, v in point.items() if k != "lma"}
    assert model_logp(model, free) == pytest.approx(expected - prior_logp("lma", 185.0), rel=1e-9)
    with pytest.raises(ValueError, match="fixed"):
        build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT, fixed={"not_a_parameter": 1.0})


def test_cold_limit_model_equals_hand_computation():
    """D55: with the cold limit, the log posterior still equals the hand computation
    (numpy model with the cold limit, plus the two threshold priors)."""
    observed, sigma, use = fake_observations(median_parameters())
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True)
    expected, point = hand_logp(median_parameters(), 0.2, 0.4, observed, sigma, use,
                                "randunc_plus_model_error", cold=(-6.0, 5.0))
    assert model_logp(model, point) == pytest.approx(expected, rel=1e-9)


def test_fixed_parameters_with_cold_limit():
    """D56: fixing parameters and the cold limit work together: the fixed ones become
    constants, the two cold thresholds are still sampled, and the model can be evaluated."""
    observed, sigma, use = fake_observations(median_parameters())
    fixed = {"lma": 185.0, "c_lf": 0.2621, "ceff": 10.0}
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True, fixed=fixed)
    free = {rv.name for rv in model.free_RVs}
    assert not free & set(fixed)
    assert {"cold_t_shutdown_fraction", "cold_t_full_fraction"} <= free
    assert np.isfinite(model.compile_logp()(model.initial_point()))