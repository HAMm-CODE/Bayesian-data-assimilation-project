"""Is a DALEC2 run a believable boreal pine forest? The D40 checks, in one place.

Used by the prior predictive check (script 03) and by the fake-data check
(script 04, which needs a believable "true" forest and believable starting
points for the sampler, D45 and D48). One definition, so the two cannot drift
apart. Change the numbers here only together with DECISIONS.md D40.

Tier 1 = physically impossible; Tier 2 = plausible for Hyytiala.
"""

from typing import Any

import numpy as np

from dalec.model_numpy import POOLS, run
from dalec.priors import steady_state_ratio

TIER1_TARGET = 0.05   # at most 5% of prior draws physically impossible
TIER2_TARGET = 0.50   # at least half of the prior draws pass every Tier-2 check

TIER2_BANDS = {       # name: (low, high, label)
    "gpp": (500.0, 2200.0, "annual GPP (g C m$^{-2}$ yr$^{-1}$)"),
    "rh": (145.0, 740.0, "annual heterotrophic respiration (g C m$^{-2}$ yr$^{-1}$)"),
    "nee": (-700.0, 300.0, "annual NEE (g C m$^{-2}$ yr$^{-1}$)"),
    "peak_lai": (1.5, 8.0, "peak LAI (mean of yearly maxima)"),
    "foliage_ratio": (1.0 / 3.0, 3.0, "foliage pool, end / start"),
    "wood_ratio": (0.8, np.inf, "wood pool, end / start"),
    "soil_ratio": (0.7, 1.3, "soil pool, end / start"),
}


def summarise_batch(out, years: np.ndarray, weights: np.ndarray | None = None) -> dict[str, np.ndarray]:
    """Per-draw numbers needed by D40, from one batch of model runs.

    `weights` (days x 365) turns daily NEE into a mean seasonal cycle; it is
    only needed for Figure L and may be left out.
    """
    year_values = np.unique(years)
    def annual(series):   # (draws, years) of yearly totals
        return np.stack([series[:, years == y].sum(axis=1) for y in year_values], axis=1)

    gpp_year = annual(out.gpp)
    lai_year_max = np.stack([out.lai[:, years == y].max(axis=1) for y in year_values], axis=1)
    finite = np.isfinite(out.nee).all(axis=1) & np.all(
        [np.isfinite(out.pools_end[p]) for p in POOLS], axis=0)
    summary = {
        "finite": finite,
        "negative_pool": np.any([out.pools_min[p] < 0 for p in POOLS], axis=0),
        "lai_max": out.lai.max(axis=1),
        "gpp_daily_max": out.gpp.max(axis=1),
        "gpp_annual_max": gpp_year.max(axis=1),
        "gpp": gpp_year.mean(axis=1),
        "rh": annual(out.rh).mean(axis=1),
        "nee": annual(out.nee).mean(axis=1),
        "peak_lai": lai_year_max.mean(axis=1),
        "foliage_ratio": out.pools_end["c_fol"] / out.pools_start["c_fol"],
        "wood_ratio": out.pools_end["c_woo"] / out.pools_start["c_woo"],
        "soil_ratio": out.pools_end["c_som"] / out.pools_start["c_som"],
        "foliage_rs": steady_state_ratio(out.pools_in["c_fol"], out.pools_out["c_fol"],
                                         out.pools_yearly["c_fol"].mean(axis=-1),
                                         out.pools_start["c_fol"]),
    }
    if weights is not None:
        summary["nee_cycle"] = out.nee @ weights   # (draws, 365)
    return summary


def run_and_summarise(params, drivers, latitude: float, frost: float | None, years: np.ndarray,
                      weights: np.ndarray | None, batch_size: int,
                      verbose: bool = True) -> dict[str, np.ndarray]:
    """Run every draw of `params` in batches and stack the per-draw summaries."""
    p = params.as_arrays()
    n = p["f_auto"].shape[0]
    parts = []
    for start in range(0, n, batch_size):
        stop = min(start + batch_size, n)
        batch = type(params)(**{k: v[start:stop] for k, v in p.items()})
        with np.errstate(all="ignore"):   # impossible draws may overflow; they are counted
            out = run(batch, drivers, latitude_deg=latitude, frost_threshold_degc=frost)
            parts.append(summarise_batch(out, years, weights))
        if verbose:
            print(f"  draws {start + 1}-{stop} of {n} done")
    return {k: np.concatenate([part[k] for part in parts]) for k in parts[0]}


def evaluate(m: dict[str, np.ndarray]) -> dict[str, Any]:
    """True/False per draw for every D40 check."""
    checks = {
        "not a number": ~m["finite"],
        "a pool below zero": m["negative_pool"],
        "LAI above 12": m["lai_max"] > 12,
        "daily GPP above 20": m["gpp_daily_max"] > 20,
        "annual GPP above 3,000": m["gpp_annual_max"] > 3000,
    }
    tier1_fail = np.any(list(checks.values()), axis=0)
    tier2 = {name: (m[name] >= low) & (m[name] <= high) & ~tier1_fail
             for name, (low, high, _) in TIER2_BANDS.items()}
    return {"tier1_reasons": checks, "tier1_fail": tier1_fail, "tier2": tier2,
            "tier2_all": np.all(list(tier2.values()), axis=0)}