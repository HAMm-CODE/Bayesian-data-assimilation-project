#!/bin/bash
#SBATCH --job-name=dalec-fit-quick
#SBATCH --account=project_2020170
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=jobs/logs/fit-quick-%j.out
#
# Quick real-data fit (DECISIONS.md D51): 2 years, 2 short chains.
# Shows that every step runs on the real NEE. Judges nothing.
# Submit from the repository root:   sbatch jobs/fit_quick.sbatch

set -euo pipefail
source jobs/roihu_setup.sh
python scripts/05_fit_real_data.py --quick "$@"
