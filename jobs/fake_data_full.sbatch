#!/bin/bash
#SBATCH --job-name=dalec-fake-full
#SBATCH --account=project_2020170
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=jobs/logs/fake-full-%j.out
#
# Full fake-data check (DECISIONS.md D48, D49): 14 years, 4 chains of
# 1,000 warm-up + 1,000 kept draws, one chain per CPU (config sampler.cores = 4).
# Time limit 48 h: the first full run took 29.4 h (DECISIONS.md F16, D50).
# small allows up to 72 h, and only the time actually used is billed.
# Submit from the repository root:   sbatch jobs/fake_data_full.sbatch
# Speed test with a dense mass matrix (D52):
#   sbatch jobs/fake_data_full.sbatch --mass-matrix dense
# or, to start only if the quick job succeeds:
#   sbatch --dependency=afterok:<quick job id> jobs/fake_data_full.sbatch

set -euo pipefail
source jobs/roihu_setup.sh
python scripts/04_fake_data_check.py "$@"
