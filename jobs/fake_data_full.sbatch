#!/bin/bash
#SBATCH --job-name=dalec-fake-full
#SBATCH --account=project_2020170
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=jobs/logs/fake-full-%j.out
#
# Full fake-data check (DECISIONS.md D48, D49): 14 years, 4 chains of
# 1,000 warm-up + 1,000 kept draws, one chain per CPU (config sampler.cores = 4).
# Time limit 24 h: the previous repository's similar model took 3.9 h on Roihu;
# the cloud test machine projected up to about 14 h. small allows up to 72 h,
# and only the time actually used is billed.
# Submit from the repository root:   sbatch jobs/fake_data_full.sbatch
# or, to start only if the quick job succeeds:
#   sbatch --dependency=afterok:<quick job id> jobs/fake_data_full.sbatch

set -euo pipefail
source jobs/roihu_setup.sh
python scripts/04_fake_data_check.py
