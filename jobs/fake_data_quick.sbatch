#!/bin/bash
#SBATCH --job-name=dalec-fake-quick
#SBATCH --account=project_2020170
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=jobs/logs/fake-quick-%j.out
#
# Quick fake-data check (DECISIONS.md D48): 2 years, 2 short chains.
# Shows that every step runs on Roihu and prints a time estimate. Judges nothing.
# Submit from the repository root:   sbatch jobs/fake_data_quick.sbatch

set -euo pipefail
source jobs/roihu_setup.sh
python scripts/04_fake_data_check.py --quick
