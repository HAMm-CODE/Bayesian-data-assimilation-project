"""Prior sensitivity (D69): the fit name gives back the model, and the log prior
and log likelihood added to a fit are the right ones. Run with:  pytest"""

import importlib.util
from pathlib import Path

import arviz as az
import numpy as np
import pymc as pm
import pytest

from dalec.config import load_config
from dalec.inference import build_model
from dalec.priors import fixed_values, priors_for
from test_inference import DRIVERS, LAT, fake_observations
from test_model_numpy import median_parameters

_path = Path(__file__).resolve().parents[1] / "scripts" / "11_prior_sensitivity.py"
_spec = importlib.util.spec_from_file_location("prior_sensitivity_script", _path)
script = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(script)


@pytest.mark.parametrize("label, expected", [
    ("fit_model_error_MA_cold_dense_ta995_d2000",
     dict(observation_model="randunc_plus_model_error", cold=True, fix_set="MA", prior_set=None, quick=False)),
    ("fit_model_error_cold_dense_ta99",
     dict(observation_model="randunc_plus_model_error", cold=True, fix_set=None, prior_set=None, quick=False)),
    ("fit_model_error_M0_cold_lit_wood_soil_dense_ta99_quick",
     dict(cold=True, fix_set="M0", prior_set="lit_wood_soil", quick=True)),
    ("fit_cold", dict(observation_model="randunc_only", cold=True, fix_set=None)),
])
def test_options_from_label(label, expected):
    opt = script.options_from_label(label, load_config())
    for key, value in expected.items():
        assert opt[key] == value, key


def test_log_prior_and_likelihood_match_the_model():
    """On a short block: the added groups sum to the model's prior and likelihood terms."""
    priors = priors_for(cold=True)
    fixed = fixed_values(load_config()["fix_sets"]["M0"], priors)
    observed, sigma, use = fake_observations(median_parameters())
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True, fixed=fixed)
    names = [rv.name for rv in model.free_RVs]
    with model:
        draws = pm.sample_prior_predictive(draws=12, random_seed=3, var_names=names)
    idata = az.from_dict({"posterior": {n: np.asarray(draws.prior[n].values) for n in names}})
    idata, table = script.sensitivity(idata, model)
    free = [rv.name for rv in model.free_RVs]
    assert set(free) <= set(table.index)
    assert {"prior", "likelihood", "diagnosis"} <= set(table.columns)

    # one draw by hand: the log prior summed over quantities and the log likelihood
    # summed over days equal the model's own terms at that point
    point = {}
    for rv in model.free_RVs:
        value = float(np.asarray(idata.posterior[rv.name].values).reshape(-1)[0])
        point[rv.name] = value
    ip = model.initial_point()
    for rv in model.free_RVs:
        transform = model.rvs_to_transforms[rv]
        v = point[rv.name]
        ip[model.rvs_to_values[rv].name] = (transform.forward(v, *rv.owner.inputs).eval() if transform else v)
    prior_terms = model.compile_fn(pm.math.sum([pm.math.sum(t) for t in model.logp(vars=model.free_RVs, jacobian=False, sum=False)]),
                                   on_unused_input="ignore")(ip)
    like_terms = model.compile_fn(model.logp(vars=model.observed_RVs), on_unused_input="ignore")(ip)
    got_prior = sum(float(np.asarray(idata.log_prior[n].values).reshape(-1)[0]) for n in free)
    got_like = float(np.asarray(idata.log_likelihood["nee"].values).reshape(idata.log_likelihood["nee"].shape[0] * idata.log_likelihood["nee"].shape[1], -1)[0].sum())
    assert got_prior == pytest.approx(float(prior_terms), rel=1e-8)
    assert got_like == pytest.approx(float(like_terms), rel=1e-8)


def test_own_prior_table_has_one_row_per_quantity():
    priors = priors_for(cold=True)
    fixed = fixed_values(load_config()["fix_sets"]["M0"], priors)
    observed, sigma, use = fake_observations(median_parameters())
    model = build_model(DRIVERS, observed, sigma, use, latitude_deg=LAT,
                        observation_model="randunc_plus_model_error", cold_limit=True, fixed=fixed)
    names = [rv.name for rv in model.free_RVs]
    with model:
        draws = pm.sample_prior_predictive(draws=40, random_seed=4, var_names=names)
    idata = az.from_dict({"posterior": {n: np.asarray(draws.prior[n].values) for n in names}})
    idata, _ = script.sensitivity(idata, model)
    own = script.own_prior_table(idata, names)
    assert list(own.index) == names
    assert {"own_prior", "likelihood", "reading"} <= set(own.columns)
    assert np.all(np.isfinite(own[["own_prior", "likelihood"]].to_numpy(dtype=float)))
