"""Turning numbers from another DALEC2 analysis into priors (DECISIONS.md D64).

The source is the CARDAMOM 2001-2010 analysis of Bloom et al. (2016), read at
the 1 degree grid cell of FI-Hyy. The rules were fixed in D64 before any value
was read:
    centre  = the cell's median;
    spread  = the cell's spread on the same scale, doubled (the other analysis
              counts with a quarter of its weight: a "power prior" with weight 1/4);
    rates   = 1 / (365.25 x residence time in years), and for the soil also
              divided by the site's mean temperature factor, because our soil
              rate is the rate at 0 degC (D20, D35).
Each function returns a prior in the format of `priors.PRIORS`.
"""

import numpy as np

WIDEN = 2.0            # D64: spreads doubled
DAYS_PER_YEAR = 365.25


def lognormal_prior(median: float, sdlog: float, widen: float = WIDEN) -> tuple:
    """Log-normal prior: median as given, log sd widened."""
    if not (median > 0 and sdlog > 0):
        raise ValueError(f"need a positive median and log sd, got {median}, {sdlog}")
    return ("lognormal", float(median), float(widen * sdlog))


def rate_prior(residence_years: float, sdlog: float, temperature_factor: float = 1.0,
               widen: float = WIDEN) -> tuple:
    """Log-normal prior on a turnover rate (d-1) from a residence time (years).
    On the log scale a rate is minus the log residence time, so the log sd is the
    same. `temperature_factor` turns a residence time that includes the effect of
    temperature into the rate at 0 degC (D64 rule 3)."""
    if not (residence_years > 0 and temperature_factor > 0):
        raise ValueError("residence time and temperature factor must be positive")
    rate = 1.0 / (DAYS_PER_YEAR * residence_years * temperature_factor)
    return lognormal_prior(rate, sdlog, widen)


def beta_prior(mean: float, sd: float, widen: float = WIDEN) -> tuple:
    """Beta prior with the given mean and the sd widened (D64 rule 5)."""
    sd = widen * sd
    if not (0 < mean < 1) or sd <= 0 or sd**2 >= mean * (1 - mean):
        raise ValueError(f"no beta has mean {mean} and sd {sd}")
    k = mean * (1 - mean) / sd**2 - 1
    return ("beta", float(mean * k), float((1 - mean) * k))


def dirichlet_prior(shares_of_gpp: list[float], sds_of_gpp: list[float], f_auto: float,
                    widen: float = WIDEN) -> tuple:
    """Dirichlet prior on the four growth shares (labile, foliage, roots, wood),
    from their shares of GPP (D64 rule 5). Means: each share divided by the sum of
    the four. Concentration: chosen so that the average sd of the four shares of
    growth equals the average widened sd (each sd of GPP divided by 1 - f_auto).
    For a Dirichlet with means m_i and concentration c, sd_i^2 = m_i (1 - m_i) / (c + 1)."""
    shares = np.asarray(shares_of_gpp, dtype=float)
    sds = widen * np.asarray(sds_of_gpp, dtype=float) / (1.0 - f_auto)
    if shares.shape != (4,) or sds.shape != (4,) or np.any(shares <= 0) or np.any(sds <= 0):
        raise ValueError("need four positive shares and four positive sds")
    means = shares / shares.sum()
    # c + 1 = (mean of sqrt(m (1 - m)) / mean sd)^2 matches the average sd exactly
    target = sds.mean()
    c = (np.sqrt(means * (1 - means)).mean() / target) ** 2 - 1.0
    if c <= 0:
        raise ValueError("the shares are too uncertain for a Dirichlet with these means")
    return ("dirichlet", *[float(c * m) for m in means])


def sdlog_from_percentiles(p_low: float, p_high: float, z: float) -> float:
    """Log sd from two percentiles that are z standard deviations either side of
    the median (for example the 25th and 75th: z = 0.6745)."""
    return float((np.log(p_high) - np.log(p_low)) / (2.0 * z))


def step_fraction_priors(weights: tuple) -> dict[str, tuple]:
    """The same allocation prior written as the three step fractions of D67:
    foliage first, then labile, then roots (wood takes the rest). `weights` is the
    output of `dirichlet_prior` (labile, foliage, roots, wood). The step form is the
    standard equivalent: each step is a beta of its own weight against the sum of
    the weights still to come."""
    _, w_lab, w_fol, w_roo, w_woo = weights
    return {"frac_fol": ("beta", float(w_fol), float(w_lab + w_roo + w_woo)),
            "frac_lab": ("beta", float(w_lab), float(w_roo + w_woo)),
            "frac_roo": ("beta", float(w_roo), float(w_woo))}
