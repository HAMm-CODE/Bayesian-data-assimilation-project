"""DALEC2 forward model in plain numpy: the readable reference version.

Follows the official DALEC2 code (CARDAMOM `DALEC_400.c`) and the decisions in
DECISIONS.md D20 to D31; the equations are written out in MODEL.md Section 3.

One day at a time, it moves carbon between six pools and returns the daily
fluxes. Every parameter may be a single number or an array with one value per
parameter draw, so 1,000 draws run together in one loop over the days. This is
the version you can read and check by hand; the fast PyTensor version for NUTS
must later reproduce it exactly.

Units: pools g C m-2, fluxes g C m-2 d-1, temperatures degC, radiation
MJ m-2 d-1, CO2 umol mol-1, rates d-1.
"""

from dataclasses import dataclass, fields
from typing import Any

import numpy as np
import pandas as pd

from dalec.acm import DALEC2_ACM, AcmConstants, frost_multiplier, gpp as acm_gpp

# W m-2 (daily mean) to MJ m-2 d-1: 86,400 s per day / 1,000,000 J per MJ.
W_M2_TO_MJ_M2_DAY = 0.0864

POOLS = ("c_lab", "c_fol", "c_roo", "c_woo", "c_lit", "c_som")

# Leaf-timing constants (MODEL.md 3.3). s = 365.25 / pi makes each pulse repeat
# every 365.25 days.
_S = 365.25 / np.pi
# Coefficients of the polynomial that gives the phase offset psi (D31), from
# DALEC_ALL/offset.c, highest power first.
_OFFSET_POLY = (0.000023599784710, 0.000332730053021, 0.000901865258885,
                -0.005437736864888, -0.020836027517787, 0.126972018064287,
                -0.188459767342504)
# The labile release uses the same formula with a fixed lifespan of 1.001
# years (official code), which empties the labile pool each spring.
_ONSET_LIFESPAN = 1.001


# ---------------------------------------------------------------------------
# Containers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Parameters:
    """All 23 DALEC2 parameters plus the derived f_woo (MODEL.md Section 6).

    Each field is a number or a 1-D array with one value per draw.
    """

    f_auto: Any
    f_lab: Any
    f_fol: Any
    f_roo: Any
    f_woo: Any
    theta_roo: Any
    theta_woo: Any
    theta_lit: Any
    theta_som: Any
    theta_min: Any
    temperature_exponent: Any
    d_onset: Any
    cr_onset: Any
    d_fall: Any
    cr_fall: Any
    c_lf: Any
    lma: Any
    ceff: Any
    c_lab_0: Any
    c_fol_0: Any
    c_roo_0: Any
    c_woo_0: Any
    c_lit_0: Any
    c_som_0: Any

    def __post_init__(self) -> None:
        total = (np.asarray(self.f_auto) + np.asarray(self.f_lab) + np.asarray(self.f_fol)
                 + np.asarray(self.f_roo) + np.asarray(self.f_woo))
        if not np.allclose(total, 1.0, atol=1e-9):
            raise ValueError("allocation fractions must add up to 1 (D28)")
        lifespan = 1.0 / np.asarray(self.c_lf, dtype=float)
        if np.any(lifespan <= 1.001) or np.any(lifespan > 8.0):
            raise ValueError("needle lifespan 1/c_lf must lie in (1.001, 8] years, "
                             "the range of the offset polynomial (D31)")

    def as_arrays(self) -> dict[str, np.ndarray]:
        """Every field as a float array of the same shape (one entry per draw)."""
        values = {f.name: np.asarray(getattr(self, f.name), dtype=float) for f in fields(self)}
        shape = np.broadcast_shapes(*(v.shape for v in values.values()))
        return {k: np.broadcast_to(v, shape).astype(float) for k, v in values.items()}


@dataclass(frozen=True)
class Drivers:
    """Daily model inputs (D19, D20, D30). All arrays have one value per day."""

    t: np.ndarray        # days since the start of the run, 1, 2, 3, ... (leaf timing)
    doy: np.ndarray      # day of year, 1 to 366 (day length)
    t_max: np.ndarray    # degC
    t_min: np.ndarray    # degC
    sw_in: np.ndarray    # MJ m-2 d-1
    co2: np.ndarray      # umol mol-1

    def __post_init__(self) -> None:
        arrays = {f.name: np.asarray(getattr(self, f.name)) for f in fields(self)}
        lengths = {k: v.shape for k, v in arrays.items()}
        if len(set(lengths.values())) != 1 or arrays["t"].ndim != 1:
            raise ValueError(f"all drivers must be 1-D and the same length, got {lengths}")
        bad = [k for k, v in arrays.items() if not np.all(np.isfinite(v))]
        if bad:
            raise ValueError(f"drivers contain missing values: {bad}; the model cannot "
                             "step through a gap")

    @property
    def t_mid(self) -> np.ndarray:
        """Respiration temperature, the midpoint of Tmax and Tmin (D20)."""
        return 0.5 * (self.t_max + self.t_min)

    @property
    def n_days(self) -> int:
        return int(len(self.t))


def drivers_from_frame(data: pd.DataFrame, config: dict[str, Any]) -> Drivers:
    """Build Drivers from the daily table returned by the loaders in data_io.

    `data` must already contain the half-hourly extremes (attach_temperature_extremes).
    Day 1 is the first row (1 January 1997 for the calibration block; D37).
    """
    drv = config["drivers"]
    ext = drv["temperature_extremes"]
    index = pd.DatetimeIndex(data.index)
    return Drivers(
        t=np.arange(1, len(data) + 1, dtype=float),
        doy=index.dayofyear.to_numpy(dtype=float),
        t_max=data[ext["t_max"]].to_numpy(dtype=float),
        t_min=data[ext["t_min"]].to_numpy(dtype=float),
        sw_in=data[drv["model"]["sw_in"]].to_numpy(dtype=float) * W_M2_TO_MJ_M2_DAY,
        co2=data[drv["model"]["co2"]].to_numpy(dtype=float),
    )


# ---------------------------------------------------------------------------
# Leaf timing (MODEL.md 3.3; D25, D30, D31)
# ---------------------------------------------------------------------------

def psi_offset(lifespan_years):
    """Phase offset psi from the official polynomial in log(L - 1) (D31)."""
    x = np.log(np.asarray(lifespan_years, dtype=float) - 1.0)
    return np.polyval(_OFFSET_POLY, x)


def _pulse(t, anchor_day, period_days, lifespan_years):
    """Fraction of a pool released on day t, as in DALEC_400.c.

    The bell-shaped pulse repeats every 365.25 days. Its size is set by the
    lifespan L: over one year the released fractions add up to log(L / (L - 1)),
    so the pool loses about 1/L of its carbon per year. Its position is set by
    the anchor day and the offset psi.
    """
    width = period_days * np.sqrt(2.0) / 2.0
    size = (np.log(lifespan_years) - np.log(lifespan_years - 1.0)) / 2.0
    offset = psi_offset(lifespan_years) * width
    shape = np.sin((t - anchor_day + offset) / _S) * _S / width
    return (2.0 / np.sqrt(np.pi)) * (size / width) * np.exp(-shape**2)


def phi_onset(t, d_onset, cr_onset):
    """Share of the labile pool moved to the foliage on day t (A7)."""
    return _pulse(t, d_onset, cr_onset, _ONSET_LIFESPAN)


def phi_fall(t, d_fall, cr_fall, c_lf):
    """Share of the foliage dropped as litter on day t (A8, anchored on d_fall: D25)."""
    return _pulse(t, d_fall, cr_fall, 1.0 / np.asarray(c_lf, dtype=float))


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Output:
    """Daily results, each of shape (draws, days); pools at start and end."""

    gpp: np.ndarray
    ra: np.ndarray       # respiration by the trees (autotrophic)
    rh: np.ndarray       # respiration from litter and soil (heterotrophic)
    nee: np.ndarray      # Ra + Rh - GPP; negative = forest takes up carbon
    lai: np.ndarray
    pools_start: dict[str, np.ndarray]
    pools_end: dict[str, np.ndarray]
    total_carbon: np.ndarray   # sum of all pools, (draws, days + 1)

    @property
    def reco(self) -> np.ndarray:
        return self.ra + self.rh


def run(params: Parameters, drivers: Drivers, *, latitude_deg: float,
        frost_threshold_degc: float | None = None,
        constants: AcmConstants = DALEC2_ACM) -> Output:
    """Run DALEC2 over every day of `drivers`, for every draw in `params`.

    frost_threshold_degc: None (the main model, D24) or a temperature below
    which GPP is set to zero (the optional frost switch).
    """
    p = params.as_arrays()
    shape = p["f_auto"].shape
    n_days = drivers.n_days
    t_mid = drivers.t_mid
    frost = (np.ones(n_days) if frost_threshold_degc is None
             else frost_multiplier(drivers.t_max, drivers.t_min, frost_threshold_degc))

    pools = {name: p[f"{name}_0"].copy() for name in POOLS}
    start = {name: v.copy() for name, v in pools.items()}

    out = {k: np.empty(shape + (n_days,)) for k in ("gpp", "ra", "rh", "nee", "lai")}
    total = np.empty(shape + (n_days + 1,))
    total[..., 0] = sum(pools.values())

    for i in range(n_days):
        c_lab, c_fol, c_roo, c_woo, c_lit, c_som = (pools[n] for n in POOLS)

        # Photosynthesis from today's canopy (A12, ACM).
        lai = c_fol / p["lma"]
        gpp = acm_gpp(lai=lai, t_max=drivers.t_max[i], t_min=drivers.t_min[i],
                      sw_in=drivers.sw_in[i], co2=drivers.co2[i], doy=drivers.doy[i],
                      ceff=p["ceff"], latitude_deg=latitude_deg, constants=constants)
        gpp = gpp * frost[i]

        # Decomposition speeds up with temperature; equals 1 at 0 degC (D20).
        temp_rate = np.exp(p["temperature_exponent"] * t_mid[i])

        # Leaf timing for today (D25, D30, D31).
        onset = phi_onset(drivers.t[i], p["d_onset"], p["cr_onset"])
        fall = phi_fall(drivers.t[i], p["d_fall"], p["cr_fall"], p["c_lf"])

        # All fluxes from today's pools, then all pools updated together.
        ra = p["f_auto"] * gpp
        labile_release = c_lab * onset
        leaf_litter = c_fol * fall
        root_litter = c_roo * p["theta_roo"]
        wood_litter = c_woo * p["theta_woo"]
        rh_litter = c_lit * p["theta_lit"] * temp_rate
        rh_soil = c_som * p["theta_som"] * temp_rate
        litter_to_soil = c_lit * p["theta_min"] * temp_rate

        pools["c_lab"] = c_lab + p["f_lab"] * gpp - labile_release
        pools["c_fol"] = c_fol + p["f_fol"] * gpp + labile_release - leaf_litter
        pools["c_roo"] = c_roo + p["f_roo"] * gpp - root_litter
        pools["c_woo"] = c_woo + p["f_woo"] * gpp - wood_litter
        pools["c_lit"] = c_lit + leaf_litter + root_litter - rh_litter - litter_to_soil
        pools["c_som"] = c_som + litter_to_soil + wood_litter - rh_soil

        rh = rh_litter + rh_soil
        out["gpp"][..., i] = gpp
        out["ra"][..., i] = ra
        out["rh"][..., i] = rh
        out["nee"][..., i] = ra + rh - gpp
        out["lai"][..., i] = lai
        total[..., i + 1] = sum(pools.values())

    return Output(gpp=out["gpp"], ra=out["ra"], rh=out["rh"], nee=out["nee"],
                  lai=out["lai"], pools_start=start, pools_end=pools,
                  total_carbon=total)