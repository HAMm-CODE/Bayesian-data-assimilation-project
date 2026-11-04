"""Does the prior or the data decide each posterior? Power-scaling sensitivity (DECISIONS.md D69).

For a saved real-data fit, the prior and the likelihood are each made slightly
stronger and slightly weaker (raised to a power alpha near 1), and the existing
posterior draws are re-weighted to see how much each quantity's posterior moves
(Kallioinen et al. 2024, Statistics and Computing; Gelman et al. 2020, Section 6.3).
Nothing is sampled again: the per-draw log prior and log likelihood are computed
for the saved draws, then ArviZ `psense_summary` does the re-weighting.

Reading the result, per quantity:
    prior sensitivity high, likelihood low   the prior decides (the data say little)
    both high                                the prior and the data disagree (prior-data conflict)
    both low                                 the data decide; the prior hardly matters
"High" means above 0.05, the ArviZ default.

What is scaled (D69): the prior of every sampled quantity, sigma_m included.
The soft foliage constraint (D44) is NOT scaled: it is not a prior on any one
quantity, so it is kept as it is in every perturbation.

Run from the repository root with the fit name (results/<name>.nc). The model is
rebuilt from the tags in the name, exactly as script 05 built it:
    python scripts/11_prior_sensitivity.py fit_model_error_MA_cold_dense_ta995_d2000 [--thin 2]

Saves:
    results/<fit>_psense.nc                         the fit plus its log_prior and log_likelihood groups
    figures/<fit>_prior_sensitivity.csv             one row per quantity: prior, likelihood, diagnosis
"""

import argparse
import time

import arviz as az
import numpy as np
import pymc as pm

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.inference import build_model
from dalec.model_numpy import drivers_from_frame
from dalec.plotting import heading
from dalec.priors import PRIOR_SETS, fixed_values, priors_for

_stats = getattr(pm, "stats", pm)
compute_log_prior = getattr(_stats, "compute_log_prior", None) or pm.compute_log_prior
compute_log_likelihood = getattr(_stats, "compute_log_likelihood", None) or pm.compute_log_likelihood


def options_from_label(label: str, config) -> dict:
    """The model options script 05 encoded in the fit name."""
    if "_fix-" in label:
        raise SystemExit("fits made with --fix are not supported; use a named --fix-set fit")
    tokens = label.split("_")
    fix_sets = [name for name in config["fix_sets"] if name in tokens]
    prior_sets = sorted((p for p in PRIOR_SETS if p != "site" and f"_{p}_" in f"{label}_"), key=len)
    return {
        "observation_model": ("randunc_plus_model_error" if label.startswith("fit_model_error")
                              else "randunc_only"),
        "cold": "cold" in tokens,
        "frost": "frost" in tokens,
        "free_a7": "a7" in tokens,
        "quick": tokens[-1] == "quick",
        "fix_set": fix_sets[0] if fix_sets else None,
        "prior_set": prior_sets[-1] if prior_sets else None,     # longest match: lit_wood_soil before lit
    }


def rebuild(label: str, config):
    """The same PyMC model script 05 sampled for this fit."""
    opt = options_from_label(label, config)
    if not label.startswith("fit"):
        raise SystemExit("this check is for real-data fits (names starting with 'fit')")
    priors = priors_for(opt["cold"], opt["prior_set"])
    fixed = fixed_values(config["fix_sets"][opt["fix_set"]], priors) if opt["fix_set"] else {}
    data = load_model_block(config, block="calibration")
    if opt["quick"]:
        years = data.index.year
        data = data[years < years[0] + int(config["sampler"]["quick"]["years"])]
    frost_on = opt["frost"] or bool(config["model"]["frost_cutoff"]["enabled"])
    model = build_model(
        drivers_from_frame(data, config),
        data[config["columns"]["nee"]].to_numpy(dtype=float),
        data[config["columns"]["sigma"]].to_numpy(dtype=float),
        data[MASK_USE].to_numpy(dtype=bool),
        latitude_deg=float(config["site"]["latitude_deg"]),
        frost_threshold_degc=float(config["model"]["frost_cutoff"]["threshold_degc"]) if frost_on else None,
        observation_model=opt["observation_model"],
        model_error_prior_sd=float(config["likelihood"]["model_error_prior_sd"]),
        fixed=fixed, cold_limit=opt["cold"], free_a7=opt["free_a7"], prior_set=opt["prior_set"])
    return model, opt


def sensitivity(idata, model, thin: int = 1):
    """Add log_prior and log_likelihood to the (thinned) fit and return (idata, table)."""
    if thin > 1:
        idata = idata.isel(draw=slice(None, None, thin))
    free = [rv.name for rv in model.free_RVs]
    with model:
        idata = compute_log_prior(idata, var_names=free, extend_inferencedata=True, progressbar=False)
        idata = compute_log_likelihood(idata, var_names=["nee"], extend_inferencedata=True, progressbar=False)
    table = az.psense_summary(idata, var_names=free)
    return idata, table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fit", help="real-data fit in results/, without .nc")
    parser.add_argument("--thin", type=int, default=1, help="keep every n-th draw (faster; default all)")
    args = parser.parse_args()
    config = load_config()

    heading(f"Prior sensitivity (D69): {args.fit}")
    idata = az.from_netcdf(resolve_path(config["paths"]["results"]) / f"{args.fit}.nc")
    model, opt = rebuild(args.fit, config)
    print("  model rebuilt from the name: " + ", ".join(f"{k} {v}" for k, v in opt.items()))
    began = time.perf_counter()
    idata, table = sensitivity(idata, model, args.thin)
    print(f"  log prior and log likelihood computed in {(time.perf_counter() - began) / 60:.1f} min")

    idata.to_netcdf(resolve_path(config["paths"]["results"]) / f"{args.fit}_psense.nc")
    out = resolve_path(config["paths"]["figures"]) / f"{args.fit}_prior_sensitivity.csv"
    table.to_csv(out)
    print(table.to_string())
    print(f"\n  saved {out}")
    print("  prior high + likelihood low = the prior decides; both high = prior-data conflict; "
          "both low = the data decide (threshold 0.05)")


if __name__ == "__main__":
    main()
