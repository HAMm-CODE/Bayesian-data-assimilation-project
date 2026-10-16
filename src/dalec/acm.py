"""Photosynthesis in DALEC2: the Aggregated Canopy Model (ACM).

Gives daily gross primary production (GPP, g C m-2 d-1) from the weather and
the size of the canopy. It is a line-by-line translation of the official DALEC2
code, `C/projects/CARDAMOM_MODELS/DALEC/DALEC_ALL/ACM.c` in the CARDAMOM
repository, written for numpy arrays (DECISIONS.md D21 to D23; MODEL.md 3.4).

All inputs may be numbers or numpy arrays of any shape that broadcast together
(for example: draws x days), so one call can evaluate many parameter draws at once.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AcmConstants:
    """The fixed ACM coefficients of DALEC2 (D22).

    Copied from `DALEC_400.c` (`constants[]` and `gpppars[]`). They are the same
    numbers in every DALEC2 run; they equal the "Loobos" set of Chuter et al.
    (2015, Appendix A). Not fitted here.
    """

    a2: float = 0.0156935      # day-length slope
    a3: float = 4.22273        # CO2 compensation point (umol mol-1)
    a4: float = 208.868        # CO2 half-saturation point (umol mol-1)
    a5: float = 0.0453194      # day-length offset
    a6: float = 0.37836        # hydraulic resistance coefficient
    a7: float = 7.19298        # maximum canopy light-use efficiency
    a8: float = 0.011136       # temperature coefficient, applied to Tmax
    a9: float = 2.1001         # LAI coefficient of light-use efficiency
    a10: float = 0.789798      # exponent on the water potential difference
    psi_d: float = -2.0        # soil-to-leaf water potential difference (MPa)
    r_tot: float = 1.0         # total hydraulic resistance


DALEC2_ACM = AcmConstants()


def solar_declination(doy):
    """Solar declination in radians (Chuter et al. 2015, B5; D23).

    `doy + 10` counts days from 21 December, as in `ACM.c`. Peaks at
    +23.4 degrees around 21 June (day 172 or 173).
    """
    doy = np.asarray(doy, dtype=float)
    return np.radians(-23.4 * np.cos(np.radians(360.0 * (doy + 10.0) / 365.0)))


def day_length_hours(doy, latitude_deg):
    """Hours of daylight (Chuter et al. 2015, B6).

    `ACM.c` sets 24 h when the sun never sets and 0 h when it never rises;
    clipping the arccos argument to [-1, 1] does exactly that.
    """
    mult = np.tan(np.radians(latitude_deg)) * np.tan(solar_declination(doy))
    return 24.0 * np.arccos(np.clip(-mult, -1.0, 1.0)) / np.pi


def gpp(*, lai, t_max, t_min, sw_in, co2, doy, ceff, latitude_deg, constants=DALEC2_ACM):
    """Daily GPP in g C m-2 d-1, following `ACM.c` step by step.

    Parameters (all may be arrays that broadcast together)
    ----------
    lai     : leaf area index, C_fol / lma (from the model state)
    t_max   : daily highest air temperature (degC)
    t_min   : daily lowest air temperature (degC)
    sw_in   : incoming shortwave radiation, daily total (MJ m-2 d-1)
    co2     : atmospheric CO2 (umol mol-1)
    doy     : day of year, 1 to 366 (for day length only)
    ceff    : canopy efficiency (a fitted parameter)
    latitude_deg : site latitude in degrees north
    """
    c = constants
    lai = np.asarray(lai, dtype=float)
    t_max = np.asarray(t_max, dtype=float)
    t_min = np.asarray(t_min, dtype=float)
    sw_in = np.asarray(sw_in, dtype=float)
    co2 = np.asarray(co2, dtype=float)

    # (B1) Canopy conductance: how easily CO2 enters the leaves. A wide daily
    # temperature range (dry, sunny days) lowers it. Note: HALF the range.
    g_c = np.abs(c.psi_d) ** c.a10 / (c.a6 * c.r_tot + 0.5 * (t_max - t_min))

    # (B3) What the canopy could fix if CO2 were not limiting, divided by g_c.
    # Temperature enters here, on the daily maximum.
    p = lai * ceff * np.exp(c.a8 * t_max) / g_c

    # (B2) CO2 concentration inside the leaf, from a quadratic equation.
    q = c.a3 - c.a4
    c_i = 0.5 * (co2 + q - p + np.sqrt((co2 + q - p) ** 2 - 4.0 * (co2 * q - p * c.a3)))

    # Rate at which CO2 can diffuse into the leaves.
    p_d = g_c * (co2 - c_i)

    # (B4) Canopy light-use efficiency: rises with LAI, then levels off.
    e_0 = c.a7 * lai**2 / (lai**2 + c.a9)

    # Light-limited and CO2-limited rates combined (whichever is smaller
    # dominates). With no leaves both are zero; the result is then zero, not 0/0.
    light = e_0 * sw_in
    denominator = light + p_d
    safe = np.where(denominator > 0.0, denominator, 1.0)
    p_i = np.where(denominator > 0.0, light * p_d / safe, 0.0)

    # (Eq. 12) Scale by day length.
    return p_i * (c.a2 * day_length_hours(doy, latitude_deg) + c.a5)


def frost_multiplier(t_max, t_min, threshold_degc):
    """Optional frost cut-off (D24, off in the main model).

    Returns 0 on days whose midpoint temperature is below the threshold and 1
    otherwise. Depends on the weather only, never on a parameter, so it is safe
    for gradient-based sampling. Not part of DALEC2.
    """
    t_mid = 0.5 * (np.asarray(t_max, dtype=float) + np.asarray(t_min, dtype=float))
    return np.where(t_mid < threshold_degc, 0.0, 1.0)


def cold_multiplier(t_min, t_shutdown, t_full):
    """Cold limit on photosynthesis (D55; "DALEC2cold", Stettz et al. 2022).

    0 on days whose minimum temperature is at or below `t_shutdown`, 1 at or
    above `t_full`, a straight line in between. Unlike the frost switch, the two
    thresholds are parameters (sampled in the cold-limit model). Works for
    single values or one value per draw (thresholds broadcast against t_min).
    """
    t_min = np.asarray(t_min, dtype=float)
    t_shutdown = np.asarray(t_shutdown, dtype=float)
    t_full = np.asarray(t_full, dtype=float)
    return np.clip((t_min - t_shutdown) / (t_full - t_shutdown), 0.0, 1.0)