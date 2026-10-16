#!/bin/bash
#SBATCH --job-name=dalec-fit
#SBATCH --account=project_2020170
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=jobs/logs/fit-%j.out
#
# Real-data fit (DECISIONS.md D51): 1997-2010, 4 chains of 1,000 warm-up +
# 1,000 kept draws, one chain per CPU. The fake-data check of the same size
# took 29.4 h (F16), hence 48 h.
#
# Main fit, from the repository root:
#   sbatch jobs/fit_full.sbatch
# Sensitivity runs: add the script's options after the file name, e.g.
#   sbatch jobs/fit_full.sbatch --frost
#   sbatch jobs/fit_full.sbatch --observation-model randunc_plus_model_error
# Each option set writes its own result files, so several can run at once.

set -euo pipefail
source jobs/roihu_setup.sh
python scripts/05_fit_real_data.py "$@"
