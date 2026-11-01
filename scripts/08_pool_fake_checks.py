"""Judge a model's fake-data check over two true forests together (DECISIONS.md D68).

Each run must pass the computation checks on its own (D49: no divergences, R-hat
below 1.01, effective sample size at least 400). The truth criterion is judged on
both runs pooled: all free quantities of both truths count together, and the
number of misses allowed follows the D63 rule for that pooled number (for
example at most 8 of 50 for the full cold model). Nothing is sampled.

Run from the repository root with the two run names (results/<name>.nc):
    python scripts/08_pool_fake_checks.py fake_data_check_cold_dense_ta995_d2000 \
        fake_data_check_cold_truth2_dense_ta995_d2000
"""

import argparse
import importlib.util
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd

from dalec.config import load_config, resolve_path

_spec = importlib.util.spec_from_file_location(
    "fake_check_script", Path(__file__).resolve().parent / "04_fake_data_check.py")
_fake = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fake)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs=2, help="the two fake-data check names, without .nc")
    args = parser.parse_args()
    config = load_config()
    diag, fake = config["diagnostics"], config["fake_data"]
    results, figures = resolve_path(config["paths"]["results"]), resolve_path(config["paths"]["figures"])

    all_ok, pooled = True, []
    for run in args.runs:
        table = pd.read_csv(figures / f"fake_data_recovery{run.removeprefix('fake_data_check')}.csv", index_col=0)
        free = table[table["narrowing"] != "fixed"]
        idata = az.from_netcdf(results / f"{run}.nc")
        divergences = int(np.asarray(idata.sample_stats["diverging"].values).sum())
        checks = {
            "no divergences": divergences == 0,
            f"R-hat below {diag['max_rhat']}": float(free["r_hat"].max()) < float(diag["max_rhat"]),
            f"effective sample size at least {diag['min_ess']}":
                float(min(free["ess_bulk"].min(), free["ess_tail"].min())) >= float(diag["min_ess"]),
        }
        print(f"\n{run}")
        for text, ok in checks.items():
            print(f"  {'PASS' if ok else 'FAIL'}  {text}")
        print(f"  (divergences {divergences}, largest R-hat {free['r_hat'].max():.4f}, "
              f"smallest ESS {min(free['ess_bulk'].min(), free['ess_tail'].min()):.0f}; "
              f"truth inside {int(free['truth_inside'].sum())} of {len(free)}; outside: "
              f"{', '.join(free.index[~free['truth_inside'].astype(bool)]) or 'none'})")
        all_ok = all_ok and all(checks.values())
        pooled.append(free["truth_inside"].astype(bool))

    inside = pd.concat(pooled)
    k = len(inside)
    allowed = _fake.allowed_misses(k, 1.0 - float(fake["interval"]), float(fake["miss_tail"]))
    truth_ok = int(inside.sum()) >= k - allowed
    print(f"\nPooled over both truths: inside {int(inside.sum())} of {k} "
          f"(rule: at least {k - allowed}, D63 applied to {k})  ->  {'PASS' if truth_ok else 'FAIL'}")
    print(f"\n  OVERALL (D68): {'PASS' if all_ok and truth_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
