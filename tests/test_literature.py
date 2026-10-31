"""The D64 rules for turning another analysis's numbers into priors."""

import numpy as np
import pytest
from scipy import stats

from dalec.literature import (
    beta_prior, dirichlet_prior, lognormal_prior, rate_prior, sdlog_from_percentiles,
)


def test_lognormal_doubles_the_log_sd():
    assert lognormal_prior(5000.0, 0.3) == ("lognormal", 5000.0, 0.6)


def test_rate_from_residence_time_and_temperature_factor():
    kind, rate, sdlog = rate_prior(40.0, 0.25, temperature_factor=1.2)
    assert kind == "lognormal"
    assert rate == pytest.approx(1 / (365.25 * 40.0 * 1.2))
    assert sdlog == pytest.approx(0.5)


def test_beta_has_the_mean_and_the_widened_sd():
    _, a, b = beta_prior(0.45, 0.04)
    dist = stats.beta(a, b)
    assert dist.mean() == pytest.approx(0.45)
    assert dist.std() == pytest.approx(0.08)
    with pytest.raises(ValueError):
        beta_prior(0.5, 0.3)              # sd 0.6 is impossible for a fraction


def test_dirichlet_means_and_average_sd():
    shares = [0.05, 0.10, 0.08, 0.30]     # of GPP; f_auto = 0.47
    sds = [0.01, 0.02, 0.02, 0.04]
    spec = dirichlet_prior(shares, sds, f_auto=0.47)
    alpha = np.array(spec[1:])
    m = alpha / alpha.sum()
    assert np.allclose(m, np.array(shares) / sum(shares))
    sd = np.sqrt(stats.dirichlet(alpha).var())
    assert sd.mean() == pytest.approx((2 * np.array(sds) / 0.53).mean(), rel=0.15)
    expected = np.sqrt(m * (1 - m)).mean() / np.sqrt(alpha.sum() + 1)
    assert np.sqrt(m * (1 - m) / (alpha.sum() + 1)).mean() == pytest.approx(expected)
    assert expected == pytest.approx((2 * np.array(sds) / 0.53).mean())


def test_sdlog_from_quartiles():
    x = stats.lognorm(s=0.4, scale=100.0)
    assert sdlog_from_percentiles(x.ppf(0.25), x.ppf(0.75), stats.norm.ppf(0.75)) == pytest.approx(0.4)


def test_step_fractions_reproduce_the_shares():
    """D67: the step form gives the same shares as the four-way form it replaces."""
    from dalec.literature import step_fraction_priors
    spec = ("dirichlet", 5.0, 3.0, 4.0, 20.0)          # labile, foliage, roots, wood
    steps = step_fraction_priors(spec)
    rng = np.random.default_rng(1)
    n = 200_000
    fol = rng.beta(*steps["frac_fol"][1:], n)
    lab = (1 - fol) * rng.beta(*steps["frac_lab"][1:], n)
    roo = (1 - fol - lab) * rng.beta(*steps["frac_roo"][1:], n)
    direct = rng.dirichlet(spec[1:], n)
    for got, k in ((lab, 0), (fol, 1), (roo, 2), (1 - fol - lab - roo, 3)):
        assert np.percentile(got, [5, 50, 95]) == pytest.approx(np.percentile(direct[:, k], [5, 50, 95]), abs=0.004)
