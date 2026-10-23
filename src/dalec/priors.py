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

# Cold limit on photosynthesis (D55), used only in the cold-limit model, so kept
# out of PRIORS (the prior predictive check and the main model do not use it).
# Ranges from the official code of the cold-limit DALEC (CARDAMOM v2.2,
# DALEC_831, PARS_INFO_831.c, by S. G. Stettz): lower threshold -15 to 0 degC,
# upper threshold 0 to 15 degC. The ranges do not overlap, so the upper threshold
# can never fall below the lower one. The official prior is flat over each range;
# here each threshold is a Beta(2, 2) fraction of its range (smooth, most weight
# in the middle, no hard walls for NUTS; D33).
COLD_PRIORS: dict[str, tuple] = {   # name: (range low, range high, beta a, beta b), degC
    "cold_t_shutdown": (-15.0, 0.0, 2.0, 2.0),   # no photosynthesis at or below
    "cold_t_full": (0.0, 15.0, 2.0, 2.0),        # full photosynthesis at or above
}

# D58 diagnostic: the ACM's maximum canopy light-use efficiency `a7`, fixed at the
# Loobos value 7.19298 in DALEC2 (D22), estimated instead. No boreal value exists
# to centre the prior on, so it is centred on the fixed value and kept wide:
# 95% roughly 2.7 to 19 (a factor of about 2.7 either way). Used only with
# --free-a7; kept out of PRIORS like the cold thresholds.
ACM_PRIORS: dict[str, tuple] = {
    "acm_a7": ("lognormal", 7.19298, 0.5),
}


def sample_acm(n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Draw the sampled ACM coefficient(s) of the D58 diagnostic."""
    return {name: _draw(*spec, n, rng) for name, spec in ACM_PRIORS.items()}


# Priors that differ in the cold-limit model. D57 tested returning `ceff` to its
# first prior (D38: median 16.5, sdlog 0.41) and it FAILED its prior predictive
# check (F24: 18.7% of draws with daily GPP above 20), so under the rule fixed in
# advance the cold-limit model keeps the main priors and this list is empty.
# The mechanism stays so that any later change is explicit and saved with each fit.
COLD_MODEL_PRIOR_CHANGES: dict[str, tuple] = {}


def priors_for(cold: bool = False) -> dict[str, tuple]:
    """The full prior list of one model version: the main model (PRIORS) or the
    cold-limit model (PRIORS with COLD_MODEL_PRIOR_CHANGES, D57)."""
    return {**PRIORS, **COLD_MODEL_PRIOR_CHANGES} if cold else dict(PRIORS)


def priors_to_text(priors: dict[str, tuple]) -> str:
    """The prior list as text, saved with each fit so later steps use the same priors."""
    import json
    return json.dumps(priors)


def priors_from_attrs(attrs) -> dict[str, tuple]:
    """The prior list a saved fit was made with. Fits saved before D57 carry no
    list; they all used PRIORS."""
    import json
    if "priors" not in attrs:
        return dict(PRIORS)
    return {name: tuple(spec) for name, spec in json.loads(attrs["priors"]).items()}


def sample_cold(n: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Draw the two cold-limit thresholds (D55) from COLD_PRIORS, in degC."""
    return {name: low + (high - low) * rng.beta(a, b, n)
            for name, (low, high, a, b) in COLD_PRIORS.items()}


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


SHARES_KEY = "shares"   # in a fixed set (D63): the four growth shares together, at their prior mean


def shares_mean() -> np.ndarray:
    """Prior mean of the four growth shares (shares of 1 - f_auto): alpha / sum(alpha)."""
    alpha = np.array(list(DIRICHLET_ALPHA.values()), dtype=float)
    return alpha / alpha.sum()


def fixed_values(names, priors: dict[str, tuple] | None = None) -> dict[str, Any]:
    """Values for quantities held fixed (D54, D63): each at its prior median; the four
    growth shares (named "shares") together at their prior mean, because they must sum to 1."""
    priors = PRIORS if priors is None else priors
    out: dict[str, Any] = {}
    for name in names:
        if name == SHARES_KEY:
            out[name] = shares_mean()
        elif name in priors:
            out[name] = prior_median(name, priors)
        else:
            raise ValueError(f"{name!r} cannot be fixed: it has no prior (use 'shares' for the growth shares)")
    return out


def apply_fixed(raw: dict[str, np.ndarray], fixed: dict[str, Any]) -> tuple[Parameters, dict[str, np.ndarray]]:
    """Overwrite prior draws with fixed values (D63) and rebuild the model Parameters.
    The four growth shares follow f_auto (fixed or not) times the fixed share vector."""
    raw = {k: np.array(v, dtype=float, copy=True) for k, v in raw.items()}
    n = raw["f_auto"].shape[0]
    for name, value in fixed.items():
        if name != SHARES_KEY:
            raw[name] = np.full(n, float(value))
    growth = 1.0 - raw["f_auto"]
    if SHARES_KEY in fixed:
        shares = np.broadcast_to(np.asarray(fixed[SHARES_KEY], dtype=float), (n, len(DIRICHLET_ALPHA)))
    else:   # keep each draw's own shares, but on the (possibly fixed) f_auto
        old_growth = np.sum([raw[name] for name in DIRICHLET_ALPHA], axis=0)
        shares = np.stack([raw[name] / old_growth for name in DIRICHLET_ALPHA], axis=1)
    for k, name in enumerate(DIRICHLET_ALPHA):
        raw[name] = growth * shares[:, k]
    raw["theta_min"] = raw["theta_lit"] * raw["h"] / (1.0 - raw["h"])
    return Parameters(**{k: v for k, v in raw.items() if k in Parameters.__dataclass_fields__}), raw


def sample_prior(n: int, rng: np.random.Generator,
                 priors: dict[str, tuple] | None = None,
                 fixed: dict[str, Any] | None = None) -> tuple[Parameters, dict[str, np.ndarray]]:
    """Draw n parameter sets from the priors (`PRIORS`, or another list from `priors_for`).

    Returns the Parameters for the model and a dict of every sampled quantity
    (including h and the four Dirichlet shares), for diagnostics. With `fixed`
    (from `fixed_values`), those quantities are held at their fixed values (D63);
    the random draws are made first, so the other quantities do not change.
    """
    priors = PRIORS if priors is None else priors
    raw = {name: _draw(*spec, n, rng) for name, spec in priors.items()}

    # c_lf outside the polynomial's range (about 0.3% of draws): draw again
    # from the same prior until inside. This removes a tiny tail, recorded in D41.
    low, high = C_LF_RANGE
    bad = (raw["c_lf"] < low) | (raw["c_lf"] >= high)
    while bad.any():
        raw["c_lf"][bad] = _draw(*priors["c_lf"], int(bad.sum()), rng)
        bad = (raw["c_lf"] < low) | (raw["c_lf"] >= high)

    shares = rng.dirichlet(list(DIRICHLET_ALPHA.values()), n)
    growth = 1.0 - raw["f_auto"]
    for k, name in enumerate(DIRICHLET_ALPHA):
        raw[name] = growth * shares[:, k]

    # theta_min from the share h (D35): h = theta_min / (theta_lit + theta_min)
    raw["theta_min"] = raw["theta_lit"] * raw["h"] / (1.0 - raw["h"])

    if fixed:
        return apply_fixed(raw, fixed)
    model_fields: dict[str, Any] = {k: v for k, v in raw.items() if k != "h"}
    return Parameters(**model_fields), raw


def prior_median(name: str, priors: dict[str, tuple] | None = None) -> float:
    """The median of one prior: the value a parameter is held at when it is fixed (D54).
    For `c_lf`, the median of the truncated beta (D41)."""
    from scipy import stats
    kind, a, b = (PRIORS if priors is None else priors)[name]
    if kind in ("lognormal", "normal"):
        return float(a)
    dist = stats.beta(a, b)
    if name == "c_lf":
        low, high = C_LF_RANGE
        lo, hi = dist.cdf(low), dist.cdf(high)
        return float(dist.ppf(lo + 0.5 * (hi - lo)))
    return float(dist.median())


def log_prior(raw: dict[str, Any], priors: dict[str, tuple] | None = None) -> np.ndarray:
    """Log prior density of parameter sets, the same density the PyMC model uses
    (D34 to D43, without the D44 soft term). `raw` holds the sampled quantities
    as `sample_prior` returns them: every PRIORS entry plus the four allocation
    values as shares of GPP. Used to diagnose fits (which chain is believable?)."""
    from scipy import stats
    from scipy.special import gammaln

    total = 0.0
    for name, (kind, a, b) in (PRIORS if priors is None else priors).items():
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
    if name in ACM_PRIORS:           # log scale (D58)
        return np.log(values)
    if name in COLD_PRIORS:          # a fraction of its range (D55)
        low, high = COLD_PRIORS[name][:2]
        v = np.clip((values - low) / (high - low), 1e-12, 1 - 1e-12)
        return np.log(v / (1 - v))
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
