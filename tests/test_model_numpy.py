"""Tests for src/dalec/model_numpy.py. Run from the repository root with:  pytest"""

import numpy as np
import pytest

from dalec.model_numpy import (
    Drivers, Parameters, phi_fall, phi_onset, psi_offset, run,
)

LAT = 61.8474


def synthetic_drivers(n_years=3):
    """Smooth, Hyytiala-like weather: cold dark winters, mild bright summers."""
    t = np.arange(1, 365 * n_years + 1, dtype=float)
    doy = (t - 1) % 365 + 1
    season = np.sin(2 * np.pi * (doy - 110) / 365)
    t_mid = 3.5 + 12 * season
    return Drivers(t=t, doy=doy, t_max=t_mid + 4, t_min=t_mid - 4,
                   sw_in=np.clip(10 + 10 * season, 0.5, None),
                   co2=np.full(t.size, 370.0))


def median_parameters(**changes):
    """Prior medians from DECISIONS.md D34 to D39."""
    f_auto = 0.51
    shares = np.array([13.06, 13.06, 15.26, 44.27]); shares /= shares.sum()
    f_lab, f_fol, f_roo, f_woo = (1 - f_auto) * shares
    values = dict(
        f_auto=f_auto, f_lab=f_lab, f_fol=f_fol, f_roo=f_roo, f_woo=f_woo,
        theta_roo=1.78e-3, theta_woo=5.33e-5, theta_lit=1.22e-3, theta_som=7.15e-5,
        theta_min=1.22e-3 * 0.2 / 0.8, temperature_exponent=0.0336,
        d_onset=145.0, cr_onset=40.0, d_fall=265.0, cr_fall=60.0, c_lf=0.27,
        lma=185.0, ceff=16.5,
        c_lab_0=77.0, c_fol_0=590.0, c_roo_0=238.0, c_woo_0=5900.0,
        c_lit_0=450.0, c_som_0=6100.0,
    )
    values.update(changes)
    return Parameters(**values)


# --- carbon balance -----------------------------------------------------------

def test_carbon_balance_holds_every_day():
    """MODEL.md 3.2: change in total carbon = -NEE, exactly, every day."""
    out = run(median_parameters(), synthetic_drivers(), latitude_deg=LAT)
    change = np.diff(out.total_carbon, axis=-1)
    assert np.max(np.abs(change + out.nee)) < 1e-9


def test_carbon_balance_holds_for_extreme_parameters():
    params = median_parameters(theta_lit=9e-3, theta_som=9e-4, ceff=90.0,
                               temperature_exponent=0.08, c_lf=0.6)
    out = run(params, synthetic_drivers(), latitude_deg=LAT)
    assert np.max(np.abs(np.diff(out.total_carbon, axis=-1) + out.nee)) < 1e-8


# --- leaf timing (D25, D30, D31) ----------------------------------------------

def test_offset_polynomial_reproduces_the_paper_constants():
    """D31: the polynomial at L = 1.001 gives A7's 0.6245; log(1001) gives 6.9088."""
    assert -psi_offset(1.001) / np.sqrt(2) == pytest.approx(0.6245, abs=1e-4)
    assert np.log(1001) == pytest.approx(6.9088, abs=1e-4)


def test_onset_pulse_matches_published_a7():
    """The code form equals A7 written as in the paper (MODEL.md 3.3)."""
    t = np.arange(1, 731, dtype=float)
    d, cr = 145.0, 40.0
    s = 365.25 / np.pi
    a7 = (np.sqrt(2 / np.pi) * (6.9088 / cr)
          * np.exp(-(np.sin((t - d - 0.6245 * cr) / s) * np.sqrt(2) * s / cr) ** 2))
    assert np.allclose(phi_onset(t, d, cr), a7, rtol=1e-3, atol=1e-6)


def test_labile_pool_is_emptied_by_each_spring_flush():
    t = np.arange(366, 731, dtype=float)
    left = np.prod(1 - phi_onset(t, 145.0, 40.0))
    assert left < 0.002


@pytest.mark.parametrize("c_lf, d_fall, cr_fall", [(0.27, 265, 60), (0.5, 200, 30), (0.2, 265, 100)])
def test_needle_fall_loses_about_c_lf_per_year_and_peaks_on_d_fall(c_lf, d_fall, cr_fall):
    """The pool loses about c_lf of its carbon per year, and the loss rate peaks
    on d_fall, as Bloom and Williams (2015) describe the offset."""
    t = np.arange(366, 731, dtype=float)
    phi = phi_fall(t, d_fall, cr_fall, c_lf)
    assert 1 - np.prod(1 - phi) == pytest.approx(c_lf, abs=0.03)
    pool, loss = 1.0, []
    for x in phi:
        loss.append(pool * x)
        pool *= 1 - x
    assert abs((t[int(np.argmax(loss))] - 365) - d_fall) <= 1


# --- behaviour ------------------------------------------------------------------

def test_many_draws_at_once_equal_one_at_a_time():
    ceff = np.array([8.0, 16.5, 30.0])
    together = run(median_parameters(ceff=ceff), synthetic_drivers(2), latitude_deg=LAT)
    for k, value in enumerate(ceff):
        alone = run(median_parameters(ceff=value), synthetic_drivers(2), latitude_deg=LAT)
        assert np.allclose(together.nee[k], alone.nee)


def test_median_forest_is_not_absurd():
    """Prior medians on Hyytiala-like weather: no blow-up, sensible orders of size."""
    out = run(median_parameters(), synthetic_drivers(3), latitude_deg=LAT)
    annual_gpp = out.gpp.sum() / 3
    assert 300 < annual_gpp < 3000
    assert np.all(np.isfinite(out.nee))
    assert all(np.all(v > 0) for v in out.pools_end.values())


def test_frost_switch_stops_gpp_on_cold_days():
    drivers = synthetic_drivers(1)
    out = run(median_parameters(), drivers, latitude_deg=LAT, frost_threshold_degc=-2.0)
    cold = drivers.t_mid < -2.0
    assert cold.any() and np.all(out.gpp[cold] == 0.0)


def test_rejects_allocation_that_does_not_add_up():
    with pytest.raises(ValueError, match="add up to 1"):
        median_parameters(f_woo=0.9)


def test_rejects_drivers_with_gaps():
    d = synthetic_drivers(1)
    t_max = d.t_max.copy()
    t_max[10] = np.nan
    with pytest.raises(ValueError, match="missing"):
        Drivers(t=d.t, doy=d.doy, t_max=t_max, t_min=d.t_min, sw_in=d.sw_in, co2=d.co2)


def test_pool_bookkeeping_inputs_minus_outputs_equal_change():
    """For every pool: total in - total out = end - start (needed by the D44 constraint)."""
    out = run(median_parameters(), synthetic_drivers(2), latitude_deg=LAT)
    for name in out.pools_end:
        change = out.pools_end[name] - out.pools_start[name]
        assert np.allclose(out.pools_in[name] - out.pools_out[name], change, atol=1e-8)
    # yearly snapshots on day 0, 365, 730
    assert out.pools_yearly["c_fol"].shape[-1] == 3