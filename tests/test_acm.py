"""Tests for src/dalec/acm.py. Run from the repository root with:  pytest"""

import math

import numpy as np
import pytest

from dalec.acm import DALEC2_ACM, day_length_hours, frost_multiplier, gpp, solar_declination

LAT = 61.8474  # FI-Hyy


def acm_c_reference(lai, tmax, tmin, rad, co2, doy, ceff, lat):
    """Independent scalar copy of ACM.c, line by line, using the math module.

    Written separately from acm.py so that a mistake in the numpy version
    (for example a wrong bracket or broadcasting) shows up as a mismatch.
    """
    consts = [ceff, 0.0156935, 4.22273, 208.868, 0.0453194, 0.37836, 7.19298,
              0.011136, 2.1001, 0.789798]
    pi = 3.1415927  # the C code uses this rounded pi
    gc = abs(-2.0) ** consts[9] / (consts[5] * 1.0 + 0.5 * (tmax - tmin))
    pp = lai * 1.0 / gc * consts[0] * math.exp(consts[7] * tmax)
    qq = consts[2] - consts[3]
    ci = 0.5 * (co2 + qq - pp + ((co2 + qq - pp) ** 2 - 4 * (co2 * qq - pp * consts[2])) ** 0.5)
    e0 = consts[6] * lai**2 / (lai**2 + consts[8])
    dec = -23.4 * math.cos((360.0 * (doy + 10.0) / 365.0) * pi / 180.0) * pi / 180.0
    mult = math.tan(lat * pi / 180) * math.tan(dec)
    if mult >= 1:
        dayl = 24.0
    elif mult <= -1:
        dayl = 0.0
    else:
        dayl = 24.0 * math.acos(-mult) / pi
    cps = e0 * rad * gc * (co2 - ci) / (e0 * rad + gc * (co2 - ci))
    return cps * (consts[1] * dayl + consts[4])


def test_constants_match_official_code():
    """D22: full-precision constants from DALEC_400.c."""
    c = DALEC2_ACM
    assert (c.a2, c.a3, c.a4, c.a5, c.a6) == (0.0156935, 4.22273, 208.868, 0.0453194, 0.37836)
    assert (c.a7, c.a8, c.a9, c.a10) == (7.19298, 0.011136, 2.1001, 0.789798)
    assert (c.psi_d, c.r_tot) == (-2.0, 1.0)


@pytest.mark.parametrize("case", [
    # lai, tmax, tmin, rad, co2, doy, ceff
    (3.0, 22.0, 10.0, 20.0, 370.0, 180, 16.5),   # summer day
    (3.5, -5.0, -12.0, 2.0, 380.0, 15, 16.5),    # winter day
    (1.0, 12.0, 2.0, 10.0, 360.0, 120, 7.4),     # spring, thin canopy
    (5.0, 28.0, 14.0, 28.0, 390.0, 200, 37.0),   # hot day, dense canopy
])
def test_matches_line_by_line_copy_of_acm_c(case):
    lai, tmax, tmin, rad, co2, doy, ceff = case
    expected = acm_c_reference(lai, tmax, tmin, rad, co2, doy, ceff, LAT)
    got = gpp(lai=lai, t_max=tmax, t_min=tmin, sw_in=rad, co2=co2, doy=doy,
              ceff=ceff, latitude_deg=LAT)
    # Small differences only from the C code's rounded pi (3.1415927).
    assert got == pytest.approx(expected, rel=1e-6)


def test_day_length_at_hyytiala():
    """About 19 h at midsummer, about 5 h at midwinter, about 12 h at equinox."""
    assert 18.5 < day_length_hours(172, LAT) < 19.5
    assert 4.5 < day_length_hours(355, LAT) < 5.5
    assert 11.5 < day_length_hours(80, LAT) < 12.5


def test_day_length_is_12_hours_at_equator():
    days = np.arange(1, 366)
    assert np.allclose(day_length_hours(days, 0.0), 12.0)


def test_declination_peaks_at_june_solstice():
    """D23: with doy + 10, the peak falls on day 172 or 173, not day 182."""
    days = np.arange(1, 366)
    assert days[np.argmax(solar_declination(days))] in (172, 173)


def test_polar_day_and_night_give_24_and_0_hours():
    assert day_length_hours(172, 80.0) == pytest.approx(24.0)
    assert day_length_hours(355, 80.0) == pytest.approx(0.0)


def test_no_leaves_or_no_light_gives_zero_not_nan():
    common = dict(t_max=20.0, t_min=10.0, co2=370.0, doy=180, ceff=16.5, latitude_deg=LAT)
    assert gpp(lai=0.0, sw_in=20.0, **common) == 0.0
    assert gpp(lai=3.0, sw_in=0.0, **common) == 0.0


def test_gpp_rises_with_leaves_light_and_co2():
    base = dict(lai=3.0, t_max=20.0, t_min=10.0, sw_in=20.0, co2=370.0, doy=180,
                ceff=16.5, latitude_deg=LAT)
    g0 = gpp(**base)
    assert gpp(**{**base, "lai": 4.0}) > g0
    assert gpp(**{**base, "sw_in": 25.0}) > g0
    assert gpp(**{**base, "co2": 400.0}) > g0


def test_broadcasts_over_draws_and_days():
    """One call for many parameter draws (rows) and many days (columns)."""
    n_draws, n_days = 4, 10
    ceff = np.linspace(8, 30, n_draws)[:, None]
    lai = np.full((n_draws, n_days), 3.0)
    days = np.arange(150, 150 + n_days)
    out = gpp(lai=lai, t_max=np.full(n_days, 20.0), t_min=np.full(n_days, 10.0),
              sw_in=np.full(n_days, 20.0), co2=np.full(n_days, 370.0), doy=days,
              ceff=ceff, latitude_deg=LAT)
    assert out.shape == (n_draws, n_days)
    assert np.all(np.isfinite(out)) and np.all(out > 0)
    one = gpp(lai=3.0, t_max=20.0, t_min=10.0, sw_in=20.0, co2=370.0, doy=150,
              ceff=ceff[2, 0], latitude_deg=LAT)
    assert out[2, 0] == pytest.approx(one)


def test_summer_gpp_is_in_a_believable_range():
    """A typical Hyytiala July day should give a few g C m-2 d-1, not 0.01 or 100."""
    g = gpp(lai=3.5, t_max=22.0, t_min=11.0, sw_in=20.0, co2=375.0, doy=190,
            ceff=16.5, latitude_deg=LAT)
    assert 1.0 < g < 20.0


def test_frost_multiplier():
    t_max = np.array([0.0, -1.0, 5.0])
    t_min = np.array([-6.0, -4.0, -1.0])   # midpoints -3.0, -2.5, 2.0
    assert list(frost_multiplier(t_max, t_min, -2.0)) == [0.0, 0.0, 1.0]


def test_cold_multiplier_is_a_ramp_between_the_thresholds():
    from dalec.acm import cold_multiplier
    t_min = np.array([-10.0, -7.0, -0.5, 6.0, 12.0])
    assert np.allclose(cold_multiplier(t_min, -7.0, 6.0), [0.0, 0.0, 0.5, 1.0, 1.0])
