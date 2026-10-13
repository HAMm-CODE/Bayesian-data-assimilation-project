"""DALEC2 forward model in PyTensor, for gradient-based sampling (NUTS).

The same model as `model_numpy.py`, written so that PyTensor can compute the
gradient of NEE with respect to every parameter. `tests/test_model.py` checks
that both versions give the same numbers; any change here must keep that true.

How it is organised (so the slow part is small):
1. Everything that depends only on the weather is computed once, in numpy,
   before the model is built: canopy conductance, exp(a8 Tmax), day length.
2. Everything that depends on parameters and the date, but not on the carbon
   pools, is computed for all days at once, outside the daily loop: leaf
   timing (Phi_onset, Phi_fall) and the temperature multiplier exp(Theta T).
3. Only what depends on yesterday's pools runs inside the daily loop
   (`pytensor.scan`): photosynthesis (it needs the leaf area) and the pool updates.
Parameters enter the loop as explicit `non_sequences`: the previous repository
measured a 12x slower gradient when they were reached through a Python closure.
"""

from typing import Any

import numpy as np
import pytensor
import pytensor.tensor as pt

from dalec.acm import DALEC2_ACM, AcmConstants, day_length_hours, frost_multiplier
from dalec.model_numpy import _OFFSET_POLY, _ONSET_LIFESPAN, _S, POOLS, Drivers

STEP_PARAMETERS = ("lma", "ceff", "f_auto", "f_lab", "f_fol", "f_roo", "f_woo",
                   "theta_roo", "theta_woo", "theta_lit", "theta_som", "theta_min")


# ---------------------------------------------------------------------------
# Parts computed outside the daily loop
# ---------------------------------------------------------------------------

def weather_terms(drivers: Drivers, latitude_deg: float, constants: AcmConstants = DALEC2_ACM,
                  frost_threshold_degc: float | None = None) -> dict[str, np.ndarray]:
    """ACM terms that depend on the weather only (numpy, computed once)."""
    c = constants
    g_c = np.abs(c.psi_d) ** c.a10 / (c.a6 * c.r_tot + 0.5 * (drivers.t_max - drivers.t_min))
    frost = (np.ones(drivers.n_days) if frost_threshold_degc is None
             else frost_multiplier(drivers.t_max, drivers.t_min, frost_threshold_degc))
    return {
        "g_c": g_c,
        "exp_a8_tmax": np.exp(c.a8 * drivers.t_max),
        "day_length_factor": c.a2 * day_length_hours(drivers.doy, latitude_deg) + c.a5,
        "sw_in": np.asarray(drivers.sw_in, dtype=float),
        "co2": np.asarray(drivers.co2, dtype=float),
        "frost": frost,
    }


def _psi_offset(lifespan):
    """D31 polynomial in log(L - 1), as a PyTensor expression (differentiable)."""
    x = pt.log(lifespan - 1.0)
    result = pt.zeros_like(x) + _OFFSET_POLY[0]
    for coefficient in _OFFSET_POLY[1:]:
        result = result * x + coefficient
    return result


def _pulse(t, anchor_day, period_days, lifespan):
    """Leaf-timing pulse for all days at once; same formula as model_numpy._pulse."""
    width = period_days * np.sqrt(2.0) / 2.0
    size = (pt.log(lifespan) - pt.log(lifespan - 1.0)) / 2.0
    offset = _psi_offset(lifespan) * width
    shape = pt.sin((t - anchor_day + offset) / _S) * _S / width
    return (2.0 / np.sqrt(np.pi)) * (size / width) * pt.exp(-shape**2)


# ---------------------------------------------------------------------------
# The daily loop
# ---------------------------------------------------------------------------

def _step(g_c, exp_a8_tmax, dlf, sw_in, co2, frost, onset, fall, temp_rate,
          c_lab, c_fol, c_roo, c_woo, c_lit, c_som,
          lma, ceff, f_auto, f_lab, f_fol, f_roo, f_woo,
          theta_roo, theta_woo, theta_lit, theta_som, theta_min,
          a3, a4, a7, a9):
    """One day: photosynthesis from today's canopy, then all pools updated together."""
    # ACM (acm.gpp), with the weather-only parts passed in.
    lai = c_fol / lma
    p = lai * ceff * exp_a8_tmax / g_c
    q = a3 - a4
    c_i = 0.5 * (co2 + q - p + pt.sqrt((co2 + q - p) ** 2 - 4.0 * (co2 * q - p * a3)))
    p_d = g_c * (co2 - c_i)
    e_0 = a7 * lai**2 / (lai**2 + a9)
    light = e_0 * sw_in
    denominator = light + p_d
    positive = pt.gt(denominator, 0.0)
    safe = pt.switch(positive, denominator, 1.0)       # no leaves: 0, not 0/0
    gpp = pt.switch(positive, light * p_d / safe, 0.0) * dlf * frost

    # Fluxes from today's pools (model_numpy.run).
    ra = f_auto * gpp
    labile_release = c_lab * onset
    leaf_litter = c_fol * fall
    root_litter = c_roo * theta_roo
    wood_litter = c_woo * theta_woo
    rh_litter = c_lit * theta_lit * temp_rate
    rh_soil = c_som * theta_som * temp_rate
    litter_to_soil = c_lit * theta_min * temp_rate

    new = (
        c_lab + f_lab * gpp - labile_release,
        c_fol + f_fol * gpp + labile_release - leaf_litter,
        c_roo + f_roo * gpp - root_litter,
        c_woo + f_woo * gpp - wood_litter,
        c_lit + leaf_litter + root_litter - rh_litter - litter_to_soil,
        c_som + litter_to_soil + wood_litter - rh_soil,
    )
    nee = ra + rh_litter + rh_soil - gpp
    foliage_in = f_fol * gpp + labile_release
    foliage_out = leaf_litter
    return (*new, gpp, nee, foliage_in, foliage_out)


def build_forward(theta: dict[str, Any], drivers: Drivers, *, latitude_deg: float,
                  frost_threshold_degc: float | None = None,
                  constants: AcmConstants = DALEC2_ACM) -> dict[str, Any]:
    """Build the DALEC2 graph for one parameter set.

    `theta` maps every field of model_numpy.Parameters to a scalar tensor (or a
    float). Returns tensors: daily `gpp` and `nee`, the pool trajectories
    `pools` (days + 1, 6; row 0 is the start), and `foliage_rs` for D44.
    """
    th = {k: pt.as_tensor_variable(np.float64(v)) if np.isscalar(v) else pt.cast(v, "float64")
          for k, v in theta.items()}
    w = weather_terms(drivers, latitude_deg, constants, frost_threshold_degc)
    t = pt.as_tensor_variable(np.asarray(drivers.t, dtype=float))
    t_mid = pt.as_tensor_variable(np.asarray(drivers.t_mid, dtype=float))

    onset = _pulse(t, th["d_onset"], th["cr_onset"], pt.as_tensor_variable(np.float64(_ONSET_LIFESPAN)))
    fall = _pulse(t, th["d_fall"], th["cr_fall"], 1.0 / th["c_lf"])
    temp_rate = pt.exp(th["temperature_exponent"] * t_mid)

    sequences = [pt.as_tensor_variable(w[k]) for k in
                 ("g_c", "exp_a8_tmax", "day_length_factor", "sw_in", "co2", "frost")]
    sequences += [onset, fall, temp_rate]
    initial = [th[f"{name}_0"] for name in POOLS]
    acm_constants = [pt.as_tensor_variable(np.float64(getattr(constants, k)))
                     for k in ("a3", "a4", "a7", "a9")]
    non_sequences = [th[k] for k in STEP_PARAMETERS] + acm_constants

    result = pytensor.scan(
        fn=_step,
        sequences=sequences,
        outputs_info=initial + [None, None, None, None],
        non_sequences=non_sequences,
        strict=True,
        return_updates=False,
    )
    *pool_paths, gpp, nee, foliage_in, foliage_out = result

    pools = pt.stack([pt.concatenate([pt.atleast_1d(start), path])
                      for start, path in zip(initial, pool_paths)], axis=1)

    # D44: Rs = (inputs / outputs) x (mean foliage on each 1 January) / (start).
    january = np.arange(0, drivers.n_days + 1, 365)
    foliage = pools[:, POOLS.index("c_fol")]
    foliage_rs = (foliage_in.sum() / foliage_out.sum()) * foliage[january].mean() / foliage[0]

    return {"gpp": gpp, "nee": nee, "pools": pools, "foliage_rs": foliage_rs,
            "lai": foliage[:-1] / th["lma"]}