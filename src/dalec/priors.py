"""The priors of DECISIONS.md D34 to D39, in one place.

`PRIORS` is the single written-down list of every prior. `sample_prior` draws
from it with numpy, for the prior predictive check (D40); the PyMC model
(`inference.py`, D47) reads the same list, so the two can never disagree.

Prior kinds (all smooth, no hard walls; D33):
    ("lognormal", median, sdlog)  log scale: spreads by factors
    ("beta", a, b)                a fraction between 0 and 1
    ("normal", mean, sd)          a position, such as a day of the year
"""

from typing import Any

import numpy as np

from dalec.model_numpy import Parameters

# Allocation of growth carbon to labile, foliage, roots, wood (D34).
DIRICHLET_ALPHA = {"f_lab": 13.06, "f_fol": 13.06, "f_roo": 15.26, "f_woo": 44.27}

PRIORS: dict[str, tuple] = {
    # group 1: allocation (D34); the four shares come from DIRICHLET_ALPHA
    "f_auto": ("beta", 25.5, 24.5),
    # group 2: turnover rates in d-1 (D35); theta_min comes from h
    "theta_roo": ("lognormal", 1.78e-3, 0.45),
    "theta_woo": ("lognormal", 5.33e-5, 0.39),
    "theta_lit": ("lognormal", 1.22e-3, 0.41),
    "theta_som": ("lognormal", 7.15e-5, 0.24),
    "h": ("beta", 4.0, 16.0),   # share of decomposing litter that becomes soil carbon
    # group 3: temperature exponent (D36)
    "temperature_exponent": ("lognormal", 0.0336, 0.35),
    # group 4: leaf timing (D37); days on the continuous counter, day 1 = 1 Jan 1997
    "d_onset": ("normal", 145.0, 15.0),
    "cr_onset": ("lognormal", 40.0, 0.35),
    "d_fall": ("normal", 265.0, 15.0),
    "cr_fall": ("lognormal", 60.0, 0.35),
    "c_lf": ("beta", 43.84, 122.84),   # round 3 (D43): needle life 3 to 5 yr; was Beta(14.51, 39.24)
    # group 5: canopy (D38)
    "lma": ("lognormal", 185.0, 0.25),
    "ceff": ("lognormal", 10.0, 0.354),   # round 2 (D42): REFLEX range 5 to 20 as 95%; was median 16.5, sdlog 0.41
    # group 6: starting pools in g C m-2 (D39)
    "c_lab_0": ("lognormal", 77.0, 0.50),
    "c_fol_0": ("lognormal", 590.0, 0.30),
    "c_roo_0": ("lognormal", 238.0, 0.30),
    "c_woo_0": ("lognormal", 5900.0, 0.25),
    "c_lit_0": ("lognormal", 450.0, 0.50),
    "c_som_0": ("lognormal", 6100.0, 0.25),
}

# Needle lifespan 1/c_lf must stay inside the range of the offset polynomial (D31).
C_LF_RANGE = (1.0 / 8.0, 1.0 / 1.001)


def _draw(kind: str, a: float, b: float, n: int, rng: np.random.Generator) -> np.ndarray:
    if kind == "lognormal":
        return a * np.exp(b * rng.standard_normal(n))
    if kind == "beta":
        return rng.beta(a, b, n)
    if kind == "normal":
        return a + b * rng.standard_normal(n)
    raise ValueError(f"unknown prior kind {kind!r}")


def sample_prior(n: int, rng: np.random.Generator) -> tuple[Parameters, dict[str, np.ndarray]]:
    """Draw n parameter sets from the priors.

    Returns the Parameters for the model and a dict of every sampled quantity
    (including h and the four Dirichlet shares), for diagnostics.
    """
    raw = {name: _draw(*spec, n, rng) for name, spec in PRIORS.items()}

    # c_lf outside the polynomial's range (about 0.3% of draws): draw again
    # from the same prior until inside. This removes a tiny tail, recorded in D41.
    low, high = C_LF_RANGE
    bad = (raw["c_lf"] < low) | (raw["c_lf"] >= high)
    while bad.any():
        raw["c_lf"][bad] = _draw(*PRIORS["c_lf"], int(bad.sum()), rng)
        bad = (raw["c_lf"] < low) | (raw["c_lf"] >= high)

    shares = rng.dirichlet(list(DIRICHLET_ALPHA.values()), n)
    growth = 1.0 - raw["f_auto"]
    for k, name in enumerate(DIRICHLET_ALPHA):
        raw[name] = growth * shares[:, k]

    # theta_min from the share h (D35): h = theta_min / (theta_lit + theta_min)
    raw["theta_min"] = raw["theta_lit"] * raw["h"] / (1.0 - raw["h"])

    model_fields: dict[str, Any] = {k: v for k, v in raw.items() if k != "h"}
    return Parameters(**model_fields), raw


def prior_median(name: str) -> float:
    """The median of one prior: the value a parameter is held at when it is fixed (D54).
    For `c_lf`, the median of the truncated beta (D41)."""
    from scipy import stats
    kind, a, b = PRIORS[name]
    if kind in ("lognormal", "normal"):
        return float(a)
    dist = stats.beta(a, b)
    if name == "c_lf":
        low, high = C_LF_RANGE
        lo, hi = dist.cdf(low), dist.cdf(high)
        return float(dist.ppf(lo + 0.5 * (hi - lo)))
    return float(dist.median())


def log_prior(raw: dict[str, Any]) -> np.ndarray:
    """Log prior density of parameter sets, the same density the PyMC model uses
    (D34 to D43, without the D44 soft term). `raw` holds the sampled quantities
    as `sample_prior` returns them: every PRIORS entry plus the four allocation
    values as shares of GPP. Used to diagnose fits (which chain is believable?)."""
    from scipy import stats
    from scipy.special import gammaln

    total = 0.0
    for name, (kind, a, b) in PRIORS.items():
        x = np.asarray(raw[name], dtype=float)
        if kind == "lognormal":
            total = total + stats.lognorm(s=b, scale=a).logpdf(x)
        elif kind == "normal":
            total = total + stats.norm(a, b).logpdf(x)
        else:
            dist = stats.beta(a, b)
            density = dist.logpdf(x)
            if name == "c_lf":          # truncated to the polynomial's range (D41)
                low, high = C_LF_RANGE
                inside = (x >= low) & (x <= high)
                density = np.where(inside, density - np.log(dist.cdf(high) - dist.cdf(low)), -np.inf)
            total = total + density
    alpha = np.array(list(DIRICHLET_ALPHA.values()))
    growth = 1.0 - np.asarray(raw["f_auto"], dtype=float)
    shares = np.stack([np.asarray(raw[k], dtype=float) / growth for k in DIRICHLET_ALPHA])
    total = total + (np.sum((alpha[:, None] - 1) * np.log(shares.reshape(4, -1)), axis=0).reshape(growth.shape)
                     + gammaln(alpha.sum()) - gammaln(alpha).sum())
    return np.asarray(total)


def to_prior_scale(name: str, values) -> np.ndarray:
    """A sampled quantity on the scale its prior lives on (D32): log for
    log-normal priors, logit for fractions (beta priors and the four allocation
    shares), unchanged for normal priors. Used to compare spreads fairly."""
    values = np.asarray(values, dtype=float)
    kind = PRIORS[name][0] if name in PRIORS else "dirichlet"
    if kind == "lognormal":
        return np.log(values)
    if kind in ("beta", "dirichlet"):
        v = np.clip(values, 1e-12, 1 - 1e-12)
        return np.log(v / (1 - v))
    return values


# ---------------------------------------------------------------------------
# Soft steady-state constraint on the foliage pool (D44).
#
# Follows the official DALEC2 code (CARDAMOM, DALEC_400/EDC2_400.c, "EDCs 7-12"):
# for a pool, Rm = total inputs / total outputs over the run, and
# Rs = Rm x (mean pool on each 1 January) / (pool at the start). A pool in
# balance has Rs = 1. The code rejects a draw when |log Rs| > log(EQF), and also
# stores a Gaussian version, -0.5 (log Rs / log EQF)^2 ("EDCPROB"). We use that
# smooth version, so NUTS can follow its gradient, on the foliage pool only.
# EQF = 2 is the code's default ("equilibrium factor").
# ---------------------------------------------------------------------------

EDC_EQF = 2.0


def steady_state_ratio(total_in, total_out, mean_january, start):
    """Rs of EDC2_400.c: 1 when the pool is in balance; above 1 when it grows, below 1 when it shrinks."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return (np.asarray(total_in) / np.asarray(total_out)) * np.asarray(mean_january) / np.asarray(start)


def edc_log_weight(rs, eqf: float = EDC_EQF):
    """Smooth penalty -0.5 (log Rs / log EQF)^2; 0 at balance, -0.5 at a factor of EQF away.
    A run that broke down (Rs not a positive number) gets weight zero (-inf)."""
    rs = np.asarray(rs, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_w = -0.5 * (np.log(rs) / np.log(eqf)) ** 2
    return np.where(np.isfinite(log_w), log_w, -np.inf)