"""Chain by chain: where did each chain go, and why? (Gelman et al. 2020, Sections 5.1 and 5.6)

When chains disagree (R-hat far above 1), there are two different stories:
  1. a chain is STUCK in a poor region, while the believable region fits the
     data at least as well (a sampling problem); or
  2. the data PREFER an unbelievable region: no believable forest fits the NEE
     as well, so the likelihood pulls the parameters away from their priors
     (a model problem: the folk theorem, Section 5.1).
For each chain this prints how well it fits the data (log-likelihood), how
believable its parameters are (log prior), the soft foliage term (D44), and
what forest it describes. The chain with the highest log posterior is the one
the model "prefers"; whether that chain is believable decides the story.

Run from the repository root, with the label of a saved fit:
    python scripts/diagnose_chains.py fit_quick
    python scripts/diagnose_chains.py fit
Reads results/<label>.nc. Light work (a few minutes); fine on a login node.
"""

import argparse

import arviz as az
import numpy as np
import pandas as pd
from scipy import stats

from dalec.config import load_config, resolve_path
from dalec.data_io import MASK_USE, load_model_block
from dalec.inference import SHARES
from dalec.model_numpy import Parameters, drivers_from_frame, run
from dalec.plausibility import evaluate, summarise_batch
from dalec.plotting import heading
from dalec.priors import PRIORS, edc_log_weight, log_prior, steady_state_ratio

SHOW = ("ceff", "lma", "c_lf", "f_auto", "f_fol", "temperature_exponent", "c_fol_0", "c_som_0")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("label", help="the fit's label, e.g. fit_quick or fit")
    parser.add_argument("--draws-per-chain", type=int, default=25)
    args = parser.parse_args()

    config = load_config()
    label = args.label
    frost = (float(config["model"]["frost_cutoff"]["threshold_degc"])
             if "_frost" in label or config["model"]["frost_cutoff"]["enabled"] else None)
    data = load_model_block(config, block="calibration")
    if label.endswith("_quick"):
        years_all = data.index.year
        data = data[years_all < years_all[0] + int(config["sampler"]["quick"]["years"])]
    drivers = drivers_from_frame(data, config)
    years = data.index.year.to_numpy()
    latitude = float(config["site"]["latitude_deg"])
    observed = data[config["columns"]["nee"]].to_numpy(dtype=float)
    sigma = data[config["columns"]["sigma"]].to_numpy(dtype=float)
    use = data[MASK_USE].to_numpy(dtype=bool)

    idata = az.from_netcdf(resolve_path(config["paths"]["results"]) / f"{label}.nc")
    post = idata.posterior
    n_chains, n_draws = np.asarray(post["f_auto"].values).shape
    pick = np.linspace(0, n_draws - 1, min(args.draws_per_chain, n_draws)).astype(int)
    sampler_lp = np.asarray(idata.sample_stats["lp"].values)

    heading(f"Chain by chain: {label} ({drivers.n_days} days, {int(use.sum())} assimilable; "
            f"{pick.size} draws per chain)")
    rows = []
    for c in range(n_chains):
        def values(name):
            return np.asarray(post[name].values, dtype=float)[c, pick]
        params = Parameters(**{f: values(f) for f in Parameters.__dataclass_fields__})
        raw = {name: values(name) for name in PRIORS}
        raw.update({name: values(name) for name in SHARES})
        out = run(params, drivers, latitude_deg=latitude, frost_threshold_degc=frost)
        sd = np.broadcast_to(sigma[use], (pick.size, int(use.sum())))
        if "sigma_m" in post.data_vars:
            sd = np.sqrt(sd**2 + values("sigma_m")[:, None] ** 2)
        log_lik = stats.norm(out.nee[:, use], sd).logpdf(observed[use]).sum(axis=1)
        rs = steady_state_ratio(out.pools_in["c_fol"], out.pools_out["c_fol"],
                                out.pools_yearly["c_fol"].mean(axis=-1), out.pools_start["c_fol"])
        soft = edc_log_weight(rs)
        prior = log_prior(raw)
        summary = summarise_batch(out, years)
        tier1 = evaluate(summary)["tier1_fail"]
        row = {
            "chain": c,
            "log-likelihood": np.median(log_lik),
            "log prior": np.median(prior),
            "soft foliage term": np.median(soft),
            "log posterior": np.median(log_lik + prior + soft),
            "sampler lp": np.median(sampler_lp[c]),
            "impossible (Tier 1)": tier1.mean(),
            "annual GPP": np.median(summary["gpp"]),
            "peak LAI": np.median(summary["peak_lai"]),
        }
        row.update({name: np.median(raw[name]) if name in raw else np.median(values(name))
                    for name in SHOW})
        rows.append(row)
    table = pd.DataFrame(rows).set_index("chain")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(table.T.to_string(float_format=lambda v: f"{v:.4g}"))

    heading("Reading it")
    best = table["log posterior"].idxmax()
    believable = table["impossible (Tier 1)"] < 0.5
    print(f"  highest log posterior: chain {best} "
          f"({'believable' if believable[best] else 'physically impossible'} forest)")
    for c in table.index:
        if c == best:
            continue
        gap_post = table.loc[best, "log posterior"] - table.loc[c, "log posterior"]
        gap_lik = table.loc[best, "log-likelihood"] - table.loc[c, "log-likelihood"]
        print(f"  chain {c}: log posterior {gap_post:.1f} lower than chain {best}; its fit to the "
              f"data (log-likelihood) is {abs(gap_lik):.1f} {'worse' if gap_lik > 0 else 'better'}")
    if not believable[best]:
        print("  -> the data prefer an impossible forest: story 2 (a model problem, Section 5.1).")
    elif not believable.all():
        print("  -> a believable chain is best, so the impossible chain(s) look stuck: story 1 "
              "(a sampling problem).")
    else:
        print("  -> every chain describes a believable forest; compare their parameters for "
              "separate modes.")
    print("  (Differences of a few units are noise; tens to thousands are decisive.)")


if __name__ == "__main__":
    main()
