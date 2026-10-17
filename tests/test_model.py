"""The PyTensor model must reproduce the numpy model exactly, and its gradients
must be right. Run from the repository root with:  pytest"""

import numpy as np
import pytensor
import pytensor.tensor as pt
import pytest

from dalec.model import build_forward
from dalec.model_numpy import POOLS, run
from dalec.priors import steady_state_ratio
from test_model_numpy import LAT, median_parameters, synthetic_drivers

DRIVERS = synthetic_drivers(2)


def compile_for(params, frost=None, outputs=("nee", "gpp", "pools", "foliage_rs")):
    graph = build_forward(params.as_arrays(), DRIVERS, latitude_deg=LAT,
                          frost_threshold_degc=frost)
    return pytensor.function([], [graph[k] for k in outputs])()


@pytest.mark.parametrize("changes", [
    {},
    {"ceff": 20.0, "lma": 140.0, "c_lf": 0.3},
    {"temperature_exponent": 0.06, "theta_som": 1.1e-4, "d_onset": 160.0, "d_fall": 250.0},
])
def test_pytensor_matches_numpy(changes):
    params = median_parameters(**changes)
    reference = run(params, DRIVERS, latitude_deg=LAT)
    nee, gpp, pools, rs = compile_for(params)
    assert np.allclose(nee, reference.nee, rtol=1e-10, atol=1e-10)
    assert np.allclose(gpp, reference.gpp, rtol=1e-10, atol=1e-10)
    for k, name in enumerate(POOLS):
        assert pools[-1, k] == pytest.approx(float(reference.pools_end[name]), rel=1e-10)
    expected_rs = steady_state_ratio(reference.pools_in["c_fol"], reference.pools_out["c_fol"],
                                     reference.pools_yearly["c_fol"].mean(axis=-1),
                                     reference.pools_start["c_fol"])
    assert rs == pytest.approx(float(expected_rs), rel=1e-10)


def test_frost_switch_matches_numpy():
    params = median_parameters()
    reference = run(params, DRIVERS, latitude_deg=LAT, frost_threshold_degc=-2.0)
    nee, *_ = compile_for(params, frost=-2.0)
    assert np.allclose(nee, reference.nee, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("name", ["ceff", "lma", "temperature_exponent", "theta_som",
                                  "c_lf", "d_fall", "c_som_0"])
def test_gradient_matches_finite_difference(name):
    """d(total NEE)/d(parameter) from PyTensor equals a numerical slope of the numpy model."""
    base = median_parameters().as_arrays()
    symbol = pt.dscalar(name)
    theta = {**base, name: symbol}
    total = build_forward(theta, DRIVERS, latitude_deg=LAT)["nee"].sum()
    gradient = pytensor.function([symbol], pt.grad(total, symbol))(float(base[name]))

    h = 1e-6 * abs(float(base[name]))
    def total_nee(value):
        return run(median_parameters(**{name: value}), DRIVERS, latitude_deg=LAT).nee.sum()
    numeric = (total_nee(float(base[name]) + h) - total_nee(float(base[name]) - h)) / (2 * h)
    assert gradient == pytest.approx(numeric, rel=1e-4, abs=1e-8)


def test_cold_limit_matches_numpy_and_has_right_gradient():
    """D55: the PyTensor cold limit equals the numpy one, and d(total NEE)/d(threshold)
    matches a numerical slope."""
    params = median_parameters()
    cold = (-5.0, 7.0)
    reference = run(params, DRIVERS, latitude_deg=LAT, cold=cold)
    graph = build_forward(params.as_arrays(), DRIVERS, latitude_deg=LAT,
                          cold=(np.float64(cold[0]), np.float64(cold[1])))
    nee = pytensor.function([], graph["nee"])()
    assert np.allclose(nee, reference.nee, rtol=1e-10, atol=1e-10)

    t0 = pt.dscalar("t0")
    total = build_forward(params.as_arrays(), DRIVERS, latitude_deg=LAT,
                          cold=(t0, np.float64(cold[1])))["nee"].sum()
    gradient = pytensor.function([t0], pt.grad(total, t0))(cold[0])
    h = 1e-5
    def total_nee(value):
        return run(params, DRIVERS, latitude_deg=LAT, cold=(value, cold[1])).nee.sum()
    numeric = (total_nee(cold[0] + h) - total_nee(cold[0] - h)) / (2 * h)
    assert gradient == pytest.approx(numeric, rel=1e-4, abs=1e-8)


def test_sampled_a7_matches_numpy_and_has_right_gradient():
    """D58: a7 given as a value gives the same NEE in numpy and PyTensor; a7 = None
    keeps the fixed DALEC2 value; d(total NEE)/d(a7) matches a numerical slope."""
    from dalec.acm import DALEC2_ACM
    params = median_parameters()
    fixed = run(params, DRIVERS, latitude_deg=LAT).nee
    same = run(params, DRIVERS, latitude_deg=LAT, a7=DALEC2_ACM.a7).nee
    assert np.allclose(same, fixed, rtol=1e-12, atol=1e-12)

    a7 = 3.5
    reference = run(params, DRIVERS, latitude_deg=LAT, cold=(-5.0, 7.0), a7=a7)
    graph = build_forward(params.as_arrays(), DRIVERS, latitude_deg=LAT,
                          cold=(np.float64(-5.0), np.float64(7.0)), a7=np.float64(a7))
    assert np.allclose(pytensor.function([], graph["nee"])(), reference.nee, rtol=1e-10, atol=1e-10)
    assert not np.allclose(reference.nee, run(params, DRIVERS, latitude_deg=LAT, cold=(-5.0, 7.0)).nee)

    a = pt.dscalar("a7")
    total = build_forward(params.as_arrays(), DRIVERS, latitude_deg=LAT, a7=a)["nee"].sum()
    gradient = pytensor.function([a], pt.grad(total, a))(a7)
    h = 1e-5
    def total_nee(value):
        return run(params, DRIVERS, latitude_deg=LAT, a7=value).nee.sum()
    numeric = (total_nee(a7 + h) - total_nee(a7 - h)) / (2 * h)
    assert gradient == pytest.approx(numeric, rel=1e-4, abs=1e-8)