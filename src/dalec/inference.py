"""The full Bayesian model in PyMC: priors, soft constraint and likelihood (D47).

    log posterior = log priors (D34 to D43, from `priors.PRIORS`)
                  + soft foliage steady-state term (D44)
                  + log likelihood of daily NEE on assimilable days (L1, D12)

In words: the priors say what we believed before seeing the tower data; the
soft term says the needle pool should neither run away nor collapse; the
likelihood scores how close the model's daily NEE is to the measured NEE,
allowing each day an error of about its RANDUNC.

The forward model is `model.build_forward` (D46), already checked against the
numpy version. This file only joins the pieces; it adds no new equations.
"""

from typing import Any

import numpy as np
import pymc as pm
import pytensor.tensor as pt

from dalec.model import build_forward
from dalec.model_numpy import Drivers
from dalec.priors import C_LF_RANGE, COLD_PRIORS, DIRICHLET_ALPHA, EDC_EQF, PRIORS

OBSERVATION_MODELS = ("randunc_only", "randunc_plus_model_error")   # L1: main, sensitivity
SHARES = tuple(DIRICHLET_ALPHA)   # f_lab, f_fol, f_roo, f_woo


def _prior_variable(name: str, spec: tuple) -> Any:
    """One PyMC random variable from a PRIORS entry. PyMC samples log-normal
    priors on the log scale and beta priors on the logit scale (D32)."""
    kind, a, b = spec
    # Plain Python numbers that fit exactly in 32 bits (such as 145.0 or 25.5)
    # would be stored by PyTensor as 32-bit numbers, and the prior density would
    # then be computed to only about 7 digits. Making them 64-bit avoids that.
    a, b = np.float64(a), np.float64(b)
    if name == "c_lf":
        # Same rule as the numpy prior draws (D41): needle life must stay inside
        # the range of the offset polynomial (D31). Written as a truncated beta,
        # which PyMC still samples on a smooth unbounded scale (no wall for NUTS).
        low, high = C_LF_RANGE
        return pm.Truncated(name, pm.Beta.dist(alpha=a, beta=b), lower=low, upper=high)
    if kind == "lognormal":
        return pm.LogNormal(name, mu=np.log(a), sigma=b)
    if kind == "beta":
        return pm.Beta(name, alpha=a, beta=b)
    if kind == "normal":
        return pm.Normal(name, mu=a, sigma=b)
    raise ValueError(f"unknown prior kind {kind!r} for {name}")


def build_model(drivers: Drivers, observed_nee: np.ndarray, sigma: np.ndarray,
                assimilable: np.ndarray, *, latitude_deg: float,
                frost_threshold_degc: float | None = None,
                observation_model: str = "randunc_only",
                model_error_prior_sd: float = 1.0,
                edc_eqf: float = EDC_EQF,
                fixed: dict[str, float] | None = None,
                cold_limit: bool = False) -> pm.Model:
    """Build the PyMC model for one block of days.

    `observed_nee`, `sigma` (RANDUNC) and `assimilable` (the D12 mask) have one
    value per day of `drivers`. The model runs through every day; only
    assimilable days enter the likelihood (D12).

    `fixed` holds parameters to keep constant instead of sampling (D54, a
    simpler model in the sense of Gelman et al. 2020, Section 5.2). They are
    stored as constants under their own names, so every later step reads them
    exactly like sampled ones.

    `cold_limit` adds the cold limit on photosynthesis (D55, Stettz et al. 2022)
    with its two thresholds sampled from COLD_PRIORS.
    """
    fixed = dict(fixed or {})
    unknown = set(fixed) - set(PRIORS)
    if unknown:
        raise ValueError(f"only parameters with a prior can be fixed; unknown: {sorted(unknown)}")
    n = drivers.n_days
    observed_nee = np.asarray(observed_nee, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    assimilable = np.asarray(assimilable, dtype=bool)
    if not (observed_nee.shape == sigma.shape == assimilable.shape == (n,)):
        raise ValueError("observed_nee, sigma and assimilable must have one value per day")
    days = np.flatnonzero(assimilable)
    if days.size == 0:
        raise ValueError("no assimilable days in this block")
    if not np.all(np.isfinite(observed_nee[days])):
        raise ValueError("an assimilable day has no NEE value; check the D12 mask")
    if not np.all(np.isfinite(sigma[days]) & (sigma[days] > 0)):
        raise ValueError("an assimilable day has a missing or non-positive sigma; check the D12 mask")
    if observation_model not in OBSERVATION_MODELS:
        raise ValueError(f"observation_model must be one of {OBSERVATION_MODELS}, "
                         f"got {observation_model!r}")

    with pm.Model() as model:
        # Priors (D34 to D43), read from the single list in priors.py.
        rv = {name: (pm.Deterministic(name, pt.as_tensor_variable(np.float64(fixed[name])))
                     if name in fixed else _prior_variable(name, spec))
              for name, spec in PRIORS.items()}
        shares = pm.Dirichlet("shares", a=np.array(list(DIRICHLET_ALPHA.values())))

        # Parameters the model needs that follow from the sampled ones.
        theta = {name: value for name, value in rv.items() if name != "h"}
        growth = 1.0 - rv["f_auto"]
        for k, name in enumerate(SHARES):                       # D28
            theta[name] = pm.Deterministic(name, growth * shares[k])
        theta["theta_min"] = pm.Deterministic(                  # D35
            "theta_min", rv["theta_lit"] * rv["h"] / (1.0 - rv["h"]))

        cold = None
        if cold_limit:
            thresholds = []
            for name, (low, high, a, b) in COLD_PRIORS.items():
                fraction = pm.Beta(f"{name}_fraction", alpha=np.float64(a), beta=np.float64(b))
                thresholds.append(pm.Deterministic(
                    name, np.float64(low) + np.float64(high - low) * fraction))
            cold = tuple(thresholds)
        out = build_forward(theta, drivers, latitude_deg=latitude_deg,
                            frost_threshold_degc=frost_threshold_degc, cold=cold)

        # Soft steady-state constraint on the foliage pool (D44).
        rs = pm.Deterministic("foliage_rs", out["foliage_rs"])
        pm.Potential("foliage_steady_state", -0.5 * (pt.log(rs) / np.log(edc_eqf)) ** 2)

        # Likelihood (L1): NEE_t ~ Normal(model_t, sigma_t) on assimilable days.
        error_sd = pt.as_tensor_variable(sigma[days])
        if observation_model == "randunc_plus_model_error":
            sigma_m = pm.HalfNormal("sigma_m", sigma=np.float64(model_error_prior_sd))
            error_sd = pt.sqrt(error_sd**2 + sigma_m**2)
        pm.Normal("nee", mu=out["nee"][days], sigma=error_sd, observed=observed_nee[days])
    return model


def initial_values(raw: dict[str, np.ndarray], i: int) -> dict[str, Any]:
    """Starting point for one chain from prior draw number `i` of `priors.sample_prior` (D45).

    The four allocation values in `raw` are shares of GPP; the model samples them
    as shares of growth (1 - f_auto), so they are divided back here.
    """
    start = {name: float(raw[name][i]) for name in PRIORS}
    growth = 1.0 - start["f_auto"]
    start["shares"] = np.array([raw[name][i] / growth for name in SHARES])
    return start


MASS_MATRIX_INIT = {"diagonal": "adapt_diag", "dense": "adapt_full"}   # D49, D52


def sample(model: pm.Model, starts: list[dict[str, Any]], *, draws: int, tune: int,
           target_accept: float, seed: int, cores: int, mass_matrix: str = "diagonal"):
    """Run NUTS (D49): one chain per starting point, no random jitter on the start.

    PyMC normally shakes each starting point by up to a factor of e on the log
    scale ("jitter"). That could move a chain out of the plausible region the
    start was chosen from, so it is switched off.

    `mass_matrix` sets how the sampler learns the shape of the posterior during
    warm-up: "diagonal" learns one scale per parameter (the validated setting,
    F16); "dense" also learns how parameters move together (D52, speed test).
    Neither changes the model, only how fast the sampler moves through it.
    """
    if mass_matrix not in MASS_MATRIX_INIT:
        raise ValueError(f"mass_matrix must be one of {tuple(MASS_MATRIX_INIT)}, got {mass_matrix!r}")
    sampled = {rv.name for rv in model.free_RVs}
    starts = [{k: v for k, v in start.items() if k in sampled} for start in starts]
    with model:
        return pm.sample(draws=draws, tune=tune, chains=len(starts), cores=min(cores, len(starts)),
                         initvals=starts, init=MASS_MATRIX_INIT[mass_matrix],
                         target_accept=target_accept, random_seed=seed, progressbar=True)