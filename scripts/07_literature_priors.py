"""Literature priors from the CARDAMOM 2001-2010 DALEC2 analysis (DECISIONS.md D64).

Reads the 1 degree grid cell that contains FI-Hyy from the files of Bloom et al.
(2016) (doi:10.7488/ds/316, CC-BY 4.0) and turns its numbers into priors by the
rules fixed in D64 before any value was read (`src/dalec/literature.py`).
Prints the cell's numbers, the priors they give (to be written into
`PRIOR_SETS["lit"]` in `src/dalec/priors.py`) and, for comparison only, other
outputs of the same cell. Nothing is sampled; no tower data are used.

Run from the repository root, with the folder of the unpacked .nc files:
    python -I scripts/07_literature_priors.py <folder> [--temperature-factor 1.2094]

The temperature factor is the 1997-2010 mean of exp(0.0336 T_mid) at FI-Hyy
(D64 rule 3; printed on Roihu on 2026-10-27).
"""

import argparse
import sys
from pathlib import Path

import netCDF4
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))   # works with python -I

from dalec.literature import (   # noqa: E402
    beta_prior, dirichlet_prior, lognormal_prior, rate_prior, sdlog_from_percentiles,
    step_fraction_priors,
)

SITE = (61.8474, 24.2948)       # FI-Hyy latitude, longitude (degrees)
Z_QUARTILE = 0.6744897501960817  # the 25th and 75th percentiles are this many sd from the median
TC_HA_TO_G_M2 = 100.0           # 1 t C ha-1 = 100 g C m-2


def open_file(folder: Path, key: str):
    for name in (f"CARDAMOM_2001_2010_{key}.nc", f"CARDAMOM_2001_2010_{key}"):
        path = folder / name
        if path.exists():
            return netCDF4.Dataset(path)
    raise FileNotFoundError(f"no CARDAMOM file for {key} in {folder}")


def cell_index(ds) -> tuple[int, int]:
    """The grid cell whose 1 degree box contains the site (cell centres are at x.5)."""
    lat, lon = np.asarray(ds["latitude"][:]), np.asarray(ds["longitude"][:])
    i, j = int(np.argmin(np.abs(lat - SITE[0]))), int(np.argmin(np.abs(lon - SITE[1])))
    if abs(lat[i] - SITE[0]) > 0.5 or abs(lon[j] - SITE[1]) > 0.5:
        raise RuntimeError("the site is not inside the nearest cell")
    return i, j


def cell(folder: Path, key: str) -> dict[str, float]:
    """Every statistic of one file at the site's cell."""
    ds = open_file(folder, key)
    i, j = cell_index(ds)
    out = {name: float(np.ma.filled(var[i, j], np.nan)) for name, var in ds.variables.items()
           if var.dimensions == ("lat", "lon")}
    out["_units"] = next(getattr(v, "units", "") for v in ds.variables.values() if v.dimensions == ("lat", "lon"))
    out["_centre"] = (float(ds["latitude"][i]), float(ds["longitude"][j]))
    return out


def log_sd_meaning(folder: Path) -> str:
    """What the RT files' "log-based" sd is, judged over all land cells (not at the
    site): the sd of log(residence time), or exp of it (a multiplicative factor).
    Compared with the log sd implied by the quartiles of the same files."""
    ds = open_file(folder, "RT_WOO")
    sd = np.ma.filled(ds["Standard_Deviation_(log_based)"][:], np.nan)
    p25 = np.ma.filled(ds["25th_percentile"][:], np.nan).astype(float)
    p75 = np.ma.filled(ds["75th_percentile"][:], np.nan).astype(float)
    good = np.isfinite(p25) & np.isfinite(p75) & (p25 > 0) & (p75 > p25)
    q = np.full(p25.shape, np.nan)
    q[good] = (np.log(p75[good]) - np.log(p25[good])) / (2.0 * Z_QUARTILE)
    ok = good & np.isfinite(sd) & (sd > 0)
    as_log = np.median(np.abs(np.log(sd[ok] / q[ok])))
    as_factor = np.median(np.abs(np.log(np.log(np.clip(sd[ok], 1.0 + 1e-9, None)) / q[ok])))
    print(f"  (over {int(ok.sum())} land cells, median |log ratio| to the quartile log sd: "
          f"read as sd of log {as_log:.3f}, read as a factor {as_factor:.3f})")
    return "log" if as_log <= as_factor else "factor"


def rt_sdlog(values: dict[str, float], meaning: str) -> float:
    sd = values["Standard_Deviation_(log_based)"]
    return float(sd if meaning == "log" else np.log(sd))


def quartile_sdlog(values: dict[str, float]) -> float:
    return sdlog_from_percentiles(values["25th_percentile"], values["75th_percentile"], Z_QUARTILE)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("folder", type=Path, help="folder with the CARDAMOM .nc files")
    parser.add_argument("--temperature-factor", type=float, default=1.2094,
                        help="1997-2010 mean of exp(0.0336 T_mid) at FI-Hyy (D64 rule 3)")
    args = parser.parse_args()
    folder = args.folder

    meaning = log_sd_meaning(folder)
    v = {key: cell(folder, key) for key in
         ("AF_AUR", "AF_LAB", "AF_FOL", "AF_ROO", "AF_WOO", "RT_WOO", "RT_SOM", "LCMA",
          "CP_WOO", "CP_SOM")}
    print(f"cell centre {v['RT_WOO']['_centre']} (1 degree box around it contains FI-Hyy {SITE})")
    print(f"RT files: the 'log-based' sd is read as {'the sd of log' if meaning == 'log' else 'exp(sd of log)'} "
          "(judged over all land cells)")
    print("\nCARDAMOM numbers at the cell (median; spread used):")
    for key, values in v.items():
        print(f"  {key:7s} median {values['Median']:.6g} {values['_units']:12s} "
              f"5% {values['5th_percentile']:.4g}  25% {values['25th_percentile']:.4g}  "
              f"75% {values['75th_percentile']:.4g}  95% {values['95th_percentile']:.4g}  "
              + (f"log sd {rt_sdlog(values, meaning):.3f}" if key.startswith("RT")
                 else f"sd {values['Standard_Deviation']:.4g}  quartile log sd {quartile_sdlog(values):.3f}"))

    f_auto = v["AF_AUR"]["Median"]
    priors = {
        "f_auto": beta_prior(f_auto, v["AF_AUR"]["Standard_Deviation"]),
        "theta_woo": rate_prior(v["RT_WOO"]["Median"], rt_sdlog(v["RT_WOO"], meaning)),
        "theta_som": rate_prior(v["RT_SOM"]["Median"], rt_sdlog(v["RT_SOM"], meaning),
                                temperature_factor=args.temperature_factor),
        "lma": lognormal_prior(v["LCMA"]["Median"], quartile_sdlog(v["LCMA"])),
        "c_woo_0": lognormal_prior(TC_HA_TO_G_M2 * v["CP_WOO"]["Median"], quartile_sdlog(v["CP_WOO"])),
        "c_som_0": lognormal_prior(TC_HA_TO_G_M2 * v["CP_SOM"]["Median"], quartile_sdlog(v["CP_SOM"])),
        **step_fraction_priors(dirichlet_prior([v[k]["Median"] for k in ("AF_LAB", "AF_FOL", "AF_ROO", "AF_WOO")],
                                  [v[k]["Standard_Deviation"] for k in ("AF_LAB", "AF_FOL", "AF_ROO", "AF_WOO")],
                                  f_auto=f_auto)),
    }
    print("\nLiterature priors (D64), for PRIOR_SETS['lit']:")
    for name, spec in priors.items():
        print(f'    "{name}": {tuple(round(x, 6) if isinstance(x, float) else x for x in spec)},')

    print("\nFor comparison only (not used as priors):")
    for key, scale, unit in (("CP_FOL", TC_HA_TO_G_M2, "g C m-2"), ("CP_LAB", TC_HA_TO_G_M2, "g C m-2"),
                             ("CP_ROO", TC_HA_TO_G_M2, "g C m-2"), ("CP_LIT", TC_HA_TO_G_M2, "g C m-2"),
                             ("RT_FOL", 1.0, "yr"), ("RT_ROO", 1.0, "yr"), ("RT_LIT", 1.0, "yr"),
                             ("FL_GPP", 365.25, "g C m-2 yr-1"), ("FL_RHE", 365.25, "g C m-2 yr-1"),
                             ("FL_NEE", 365.25, "g C m-2 yr-1")):
        values = cell(folder, key)
        q = (f"  (25% {scale * values['25th_percentile']:.4g}, 75% {scale * values['75th_percentile']:.4g})"
             if "25th_percentile" in values else "")
        print(f"  {key:7s} median {scale * values['Median']:.4g} {unit}{q}")
    lai = v["LCMA"]["Median"] and cell(folder, "CP_FOL")["Median"] * TC_HA_TO_G_M2 / v["LCMA"]["Median"]
    print(f"  implied mean LAI (foliage / LMA, medians): {lai:.2f}")


if __name__ == "__main__":
    main()
