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

# The four shares of GPP that go to the living pools (D28, D34). They are built step
# by step, as in the official DALEC2 code (D67): foliage takes a fraction of the growth
# carbon (GPP - Ra), labile a fraction of what is left, roots a fraction of what is left
# after that, and wood the rest. Each fraction has a beta prior.
SHARE_NAMES = ("f_lab", "f_fol", "f_roo", "f_woo")
FRACTIONS = ("frac_fol", "frac_lab", "frac_roo")

PRIORS: dict[str, tuple] = {
    # group 1: allocation (D34, D67)
    "f_auto": ("beta", 25.5, 24.5),
    "frac_fol": ("beta", 13.06, 72.59),   # share of growth carbon to foliage
    "frac_lab": ("beta", 13.06, 59.53),   # share of the rest to labile
    "frac_roo": ("beta", 15.26, 44.27),   # share of the rest to fine roots; wood takes the remainder
    # group 2: turnover rates in d-1 (D35)
    "theta_roo": ("lognormal", 1.78e-3, 0.45),
    "theta_woo": ("lognormal", 5.33e-5, 0.39),
    "theta_lit": ("lognormal", 1.22e-3, 0.41),
    "theta_som": ("lognormal", 7.15e-5, 0.24),
    # litter to soil (D66): the published range of Bloom and Williams (2015, Table 1),
    # 1e-5 to 1e-2 d-1, as the 95% range of a smooth log-normal (the D35 fallback rule)
    "theta_min": ("lognormal", 3.1623e-4, 1.7623),
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


# Prior sets (D64): changes to the site priors above, by name. "site" is the
# priors as written (D34 to D39). "lit" replaces the priors of the weakly informed
# quantities by literature priors from an independent DALEC2 analysis (Bloom et
# al. 2016, CARDAMOM, the 1 degree cell of FI-Hyy). None means "not filled in yet".
PRIOR_SETS: dict[str, dict[str, tuple] | None] = {
    "site": {},
    # D64 rules applied to the CARDAMOM cell centred at 61.5 N, 24.5 E (read 2026-10-27,
    # `scripts/07_literature_priors.py`; F35). CARDAMOM medians, with the spread used:
    #   f_auto 0.498 (sd 0.092); shares of GPP labile 0.0506, foliage 0.0108, roots 0.0390,
    #   wood 0.355 (sd 0.054, 0.022, 0.061, 0.122); wood residence 30.7 yr (log sd 0.93);
    #   soil residence 54.1 yr (log sd 0.75; temperature factor 1.2094); LMA 22.8 g C m-2
    #   (quartile log sd 0.85); wood 107.6 and soil 236.1 t C ha-1 (quartile log sd 1.13, 0.39).
    # Note (F35): the step fractions have beta weights below 1 (D67 gives them in this form).
    "lit": {
        "f_auto": ("beta", 3.165244, 3.189735),
        "theta_woo": ("lognormal", 8.9e-05, 1.866693),
        "theta_som": ("lognormal", 4.2e-05, 1.499433),
        "lma": ("lognormal", 22.796788, 1.701581),
        "c_woo_0": ("lognormal", 10763.646947, 2.259382),
        "c_som_0": ("lognormal", 23613.989673, 0.775341),
        "frac_fol": ("beta", 0.006102, 0.251651),
        "frac_lab": ("beta", 0.028642, 0.223009),
        "frac_roo": ("beta", 0.02208, 0.200929),
    },
    # Exploratory, decided after reading the CARDAMOM values (F35, D65 option C): only the
    # two starting pools the supervisor named, from the same cell, so that with
    # --fix-set MA / MB they are held at the CARDAMOM medians and everything else keeps
    # the site priors.
    "lit_wood": {"c_woo_0": ("lognormal", 10763.646947, 2.259382)},
    "lit_wood_soil": {"c_woo_0": ("lognormal", 10763.646947, 2.259382),
                      "c_som_0": ("lognormal", 23613.989673, 0.775341)},
}


def priors_for(cold: bool = False, prior_set: str | None = None) -> dict[str, tuple]:
    """The full prior list of one model version: the main model (PRIORS) or the
    cold-limit model (PRIORS with COLD_MODEL_PRIOR_CHANGES, D57), optionally with
    the changes of a prior set (D64)."""
    base = {**PRIORS, **COLD_MODEL_PRIOR_CHANGES} if cold else dict(PRIORS)
    if prior_set in (None, "site"):
        return base
    if prior_set not in PRIOR_SETS:
        raise ValueError(f"unknown prior set {prior_set!r}; known: {sorted(PRIOR_SETS)}")
    changes = PRIOR_SETS[prior_set]
    if changes is None:
        raise ValueError(f"prior set {prior_set!r} is not filled in yet (D64)")
    unknown = set(changes) - set(base)
    if unknown:
        raise ValueError(f"prior set {prior_set!r} changes unknown quantities: {sorted(unknown)}")
    return {**base, **changes}


def shares_from_fractions(f_auto, frac_fol, frac_lab, frac_roo) -> dict[str, Any]:
    """The four shares of GPP (labile, foliage, roots, wood) from f_auto and the three
    step fractions (D67). Works for numbers, numpy arrays and PyTensor variables."""
    growth = 1.0 - f_auto
    f_fol = growth * frac_fol
    rest = growth - f_fol
    f_lab = rest * frac_lab
    rest = rest - f_lab
    f_roo = rest * frac_roo
    return {"f_lab": f_lab, "f_fol": f_fol, "f_roo": f_roo, "f_woo": rest - f_roo}


def fractions_from_shares(f_auto, f_lab, f_fol, f_roo) -> dict[str, Any]:
    """The step fractions back from the shares of GPP (the inverse of shares_from_fractions)."""
    growth = 1.0 - np.asarray(f_auto, dtype=float)
    after_fol = growth - np.asarray(f_fol, dtype=float)
    after_lab = after_fol - np.asarray(f_lab, dtype=float)
    return {"frac_fol": np.asarray(f_fol) / growth, "frac_lab": np.asarray(f_lab) / after_fol,
            "frac_roo": np.asarray(f_roo) / after_lab}


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


def fixed_values(names, priors: dict[str, tuple] | None = None) -> dict[str, Any]:
    """Values for quantities held fixed (D54, D63): each at its prior median."""
    priors = PRIORS if priors is None else priors
    out: dict[str, Any] = {}
    for name in names:
        if name in priors:
            out[name] = prior_median(name, priors)
        else:
            raise ValueError(f"{name!r} cannot be fixed: it has no prior")
    return out


def _model_parameters(raw: dict[str, np.ndarray]) -> Parameters:
    return Parameters(**{k: raw[k] for k in Parameters.__dataclass_fields__})


def apply_fixed(raw: dict[str, np.ndarray], fixed: dict[str, Any]) -> tuple[Parameters, dict[str, np.ndarray]]:
    """Overwrite prior draws with fixed values (D63) and rebuild the model Parameters.
    The four shares of GPP are rebuilt from f_auto and the step fractions (D67)."""
    raw = {k: np.array(v, dtype=float, copy=True) for k, v in raw.items()}
    n = raw["f_auto"].shape[0]
    for name, value in fixed.items():
        raw[name] = np.full(n, float(value))
    raw.update(shares_from_fractions(*(raw[k] for k in ("f_auto", *FRACTIONS))))
    return _model_parameters(raw), raw


def sample_prior(n: int, rng: np.random.Generator,
                 priors: dict[str, tuple] | None = None,
                 fixed: dict[str, Any] | None = None) -> tuple[Parameters, dict[str, np.ndarray]]:
    """Draw n parameter sets from the priors (`PRIORS`, or another list from `priors_for`).

    Returns the Parameters for the model and a dict of every sampled quantity
    (including the step fractions and the four shares of GPP), for diagnostics. With `fixed`
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

    raw.update(shares_from_fractions(*(raw[k] for k in ("f_auto", *FRACTIONS))))   # D67

    if fixed:
        return apply_fixed(raw, fixed)
    return _model_parameters(raw), raw


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
    as `sample_prior` returns them (every PRIORS entry, including the three step
    fractions of D67). Used to diagnose fits (which chain is believable?)."""
    from scipy import stats

    priors = PRIORS if priors is None else priors
    total = 0.0
    for name, (kind, a, b) in priors.items():
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
    return np.asarray(total)


def to_prior_scale(name: str, values) -> np.ndarray:
    """A sampled quantity on the scale its prior lives on (D32): log for
    log-normal priors, logit for fractions (beta priors and the four shares of
    GPP), unchanged for normal priors. Used to compare spreads fairly."""
    values = np.asarray(values, dtype=float)
    if name in ACM_PRIORS:           # log scale (D58)
        return np.log(values)
    if name in COLD_PRIORS:          # a fraction of its range (D55)
        low, high = COLD_PRIORS[name][:2]
        v = np.clip((values - low) / (high - low), 1e-12, 1 - 1e-12)
        return np.log(v / (1 - v))
    kind = PRIORS[name][0] if name in PRIORS else "fraction"   # the four shares of GPP
    if kind == "lognormal":
        return np.log(values)
    if kind in ("beta", "fraction"):
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
