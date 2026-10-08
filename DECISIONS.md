# Decisions

Every methodological and technical choice in this project, with its reason.
Each entry is marked with a status:

- **Decided**: settled and implemented.
- **Proposed**: preferred option chosen, not yet final.
- **Open**: still to be decided; the options are recorded.

Last updated: 2026-10-08.

---

## Part 1. Project setup

### D1. Repository content (Decided)

- Repository name: `Bayesian-data-assimilation`, named after the thesis rather than one stage of it.
- Raw and processed data are **not tracked**. The FLUXNET licence does not allow re-sharing, and everything else can be rebuilt from the raw file.
- Figures **are tracked** (small, and what reviewers look at). Sampler output (`*.nc`, `*.npz`) and cluster logs (`jobs/logs/`, `*.out`, `*.err`) are **not tracked**.
- Reason for the broad log rule: the previous repository ignored only `slurm-*.out`, so logs named `calibrate-*.out` were committed by mistake.
- `.gitignore` must be saved as UTF-8. Git cannot read a UTF-16 file (the PowerShell `>` default) and would silently ignore nothing.

### D2. Software environment (Decided)

- One dependency list only: `environment.yml`. The previous repository had four lists that disagreed with each other.
- conda with the `conda-forge` channel only (`nodefaults` blocks Anaconda's default channel), because compiled scientific libraries such as PyMC and PyTensor install most reliably this way on Windows.
- Python 3.13. Package names are unpinned while building; exact versions will be frozen into `environment.lock.yml` before producing thesis results.
- Current versions (2026-10-07): Python 3.13.16, numpy 2.5.3, pandas 3.0.6, scipy 1.18.1, matplotlib 3.11.2.
- PyMC is added when the model is built, as its own recorded change.
- On Roihu, follow CSC's current recommended way of installing Python environments; the package list stays the same.

### D7. Code layout (Decided)

- Project code lives in `src/dalec/` and is installed in editable mode (`pip install -e .`), so scripts, tests and cluster jobs import it the same way, with no path tricks.
- The `src/` layout forces Python to use the installed package, so a broken install shows up immediately.
- `pyproject.toml` holds **no dependency list**; `environment.yml` is the single list (D2). `requires-python = ">=3.13"` matches the environment actually tested.

### D8. Settings handling (Decided)

- All settings live in `config/default.yaml`. No script types a path, year or threshold directly.
- Paths resolve from the repository root (found from the location of `config.py`), not from the folder the command is run in.
- Missing or invalid settings stop the program. Nothing falls back to a guessed default (for example, there are no default years).
- The config is read with `yaml.safe_load` and UTF-8 encoding.

---

## Part 2. Data choices

### D3. Target and uncertainty columns (Decided)

- Data product: FLUXNET2015 **FULLSET daily (DD)** for FI-Hyy. SUBSET is not usable because it omits RANDUNC.
- Assimilation target: `NEE_VUT_REF`.
  - VUT (variable u* threshold): the turbulence threshold can differ between years, which suits an 18-year record better than one constant threshold (CUT).
  - REF (reference): the version most similar to all other threshold versions.
- Observation uncertainty: `NEE_VUT_REF_RANDUNC`.
- Screening: `NEE_VUT_REF_QC`, the fraction of the day that is measured or gap-filled with high confidence.
- Column names are set in the config, so switching to another version (for example `NEE_CUT_REF`) is a one-line change.

### D4. Years (Decided)

- Calibration 1997 to 2010; evaluation 2011 to 2014. Both inclusive.
- Design follows Richardson et al. (2010): one long contiguous calibration block, then a held-out block. A long record is needed because the slow pools (wood, soil) change over years to decades.
- 1996 is excluded because the CO2 driver is missing for its first 19 days, and the model cannot run through missing drivers.

### D5. QC threshold 0.75 (Decided)

- A day must have QC >= 0.75 to enter the likelihood.
- Why a threshold is needed: FLUXNET computes daily RANDUNC **from measured half-hours only**. On heavily gap-filled days, RANDUNC therefore does not describe the error of the daily value.
- Why 0.75: FLUXNET gives no recommended value. Published choices range from 0.5 to 0.8:
  - Pastorello et al. (2020): QC > 0.5, yearly data, one figure.
  - Yuan et al. (2025): QC >= 0.75, monthly data.
  - Pu et al. (2025): QC > 0.8, yearly data.
  - Tramontana et al. (2016): at most 20% low-confidence fill (equivalent to 0.8), 8-day data.
  - Kuppel et al. (2012), MDF with daily NEE: days with more than 20% missing half-hours excluded.
  - Bacour et al. (2023), MDF with daily NEE: days kept only with at least 80% of half-hours available.
- Williams et al. (2009) recommend using measured rather than gap-filled data in model-data fusion.
- Sensitivity (calibration block, 5,113 days):

  | QC threshold | Assimilable days | Share |
  |---|---|---|
  | 0.5 | 4,558 | 89.1% |
  | **0.75** | **4,413** | **86.3%** |
  | 0.9 | 4,101 | 80.2% |
  | 1.0 | 3,448 | 67.4% |

  The count is flat below 0.75 and drops steeply toward 1.0, so the choice of 0.75 changes the amount of data little. (Computed with `add_masks` at each threshold; no figure, by the author's choice.)

### D6. Seasons (Decided, 2026-10-08; currently not used)

- Chosen definition: meteorological seasons, DJF (December to February), MAM (March to May), JJA (June to August), SON (September to November).
- No current figure groups by season: the monthly box plots (Figure C) and the annual cycle (Figure E) show the seasonal pattern month by month and day by day instead. The `seasons` block was therefore removed from the config, because a setting nothing reads is misleading. Restore it if a later figure needs seasons.
- Stated limitation if used later: at Hyytiala (62 N) winter lasts longer than December to February.

### D9. Missing values and file checks (Decided)

- FLUXNET's -9999 sentinel is converted to NaN by numeric comparison (tolerance 1e-6), not by text matching, which could miss values written as `-9999.000`.
- Dates are parsed with an explicit `YYYYMMDD` format, so pandas never guesses the day/month order.
- A missing file, a missing required column, or any non-numeric column stops the program. (An earlier draft silently turned text into NaN; that was changed to follow D8.)

### D10. Calendar integrity (Decided)

- The daily file must contain every date between its first and last day exactly once, because the model steps one day per row.
- Duplicates are checked before gaps, so one missing day plus one duplicate day cannot cancel out in the row count.
- Any violation stops the program. Nothing is repaired automatically.

### D11. Year selection (Decided)

- Blocks are inclusive year ranges and must be fully covered by the file (1 January of the first year to 31 December of the last). Partial coverage stops the program.
- Calibration 1997 to 2010 = 5,113 days. Evaluation 2011 to 2014 = 1,461 days.

### D12. Assimilable days, by masking (Decided)

- No day is ever deleted. A True/False column (a mask) marks which days the likelihood may use. The forward model runs through every day, because it needs every day's drivers.
- A day is assimilable only if all three hold:
  1. NEE is present.
  2. QC >= 0.75 (D5).
  3. RANDUNC is present **and strictly positive**. A zero sigma would give infinite weight (weight = 1/sigma^2); a negative error is impossible.
- QC and RANDUNC are both needed because each removes days the other would keep (see F3).

### D13. Screening report (Decided)

- Every day is assigned to exactly one category: assimilable, or a single exclusion reason. Counts therefore never depend on the order of the filters and always sum to the total.
- The report checks itself: the rows must sum to the number of days, and the assimilable count must match the mask. A failure there means a bug in the code, not in the data.

### D14. Single entry point (Decided)

- All scripts load data through `load_screened_block()`, which reads every setting from the config.
- Lower-level functions take plain arguments and never read the config, so they can be tested in isolation.
- Sensitivity runs change a copy of the config in memory, never the file on disk.

### D15. Sigma figures (Decided, 2026-10-08)

- One script, `scripts/01_sigma_eda.py`, makes all five figures, because they share the same loaded data. Each figure is saved as its own file (PDF for the thesis, PNG for quick viewing).
- The five figures:
  - **A. Sigma over time:** RANDUNC for every day, with excluded days in grey and missing-RANDUNC periods shaded. Also prints every gap of 7 days or more.
  - **B. Sigma against |NEE|:** one colour, with a straight trend line. No fit statistics are shown, by the author's choice, to keep the figure easy to explain.
  - **C. Sigma by calendar month:** one box plot per month, all years pooled.
  - **D. Distribution of sigma:** histogram with a smoothed density curve, drawn only between the smallest and largest observed value (smoothing would otherwise suggest values below the minimum), and the median marked.
  - **E. Annual cycle of sigma:** median and middle half (25th to 75th percentile) for each day of the year. Median rather than mean, because RANDUNC is skewed. 29 February is dropped, so every year has a 365-day axis.
- Which days each figure uses (previously OPEN-2):
  - Figure A uses **every day** of the calibration block, because it shows what is excluded and missing.
  - Figures B to E use **assimilable days only**, because those are what the likelihood sees.
- Removed by the author's choice: the season weight-share figure and the QC sensitivity figure. The QC sensitivity numbers remain recorded in D5.

---

## Part 3. Likelihood

### L1. Observation model: what sigma is (Decided, 2026-10-08)

**Main model (option A): sigma is RANDUNC.**

    NEE_t ~ Normal( model_t , RANDUNC_t )      [g C m-2 d-1], assimilable days only (D12)

In words: each day's error is FLUXNET's own estimate of that day's random measurement error. Nothing is estimated on top of it.

Reasons:

- Simple and transparent: sigma comes straight from the data product, with no extra parameter to estimate or justify.
- RANDUNC is the uncertainty estimate the FLUXNET2015 pipeline provides for exactly this variable (Pastorello et al. 2020), derived following Hollinger and Richardson (2005).
- It keeps the error's dependence on flux size: days with larger fluxes get larger sigma. A constant sigma would ignore this.
- It matches the observation model of the previous repository, so no earlier result needs to be redone because of this choice.

Known consequences, stated as limitations:

- RANDUNC covers measurement random error only (computed from measured half-hours). It does not include model error, such as the known early spring onset.
- Days with very small RANDUNC (mostly winter) carry a large share of the likelihood weight (weight = 1/sigma^2; the low winter values are visible in Figures C and E), and the posterior may be narrower than the true uncertainty.
- Recent model-data fusion studies mostly use a constant sigma or one derived from model-data mismatch (Famiglietti et al. 2021; Bacour et al. 2023; Douglas et al. 2026 preprint) rather than per-day RANDUNC. Defend the choice with the reasons above and the sensitivity analysis below, not with precedent.

**Sensitivity analysis (planned, later):**

- Option C: sigma_t = sqrt( RANDUNC_t^2 + sigma_m^2 ), with sigma_m ~ HalfNormal(1) estimated jointly. Adds one shared model-error term.
- Support for option C if needed: Gelman et al. (2020, Section 10.5), golf putting example (known per-point variance plus an added independent error term "to allow for model misspecification"); Fer et al. (2018), observation-error parameters estimated during calibration.
- Compare the two on: posterior widths, winter share of the likelihood weight, R-hat, divergences and number of posterior modes.

Not chosen:

- A constant sigma (option B): ignores that sigma grows with flux size (R^2 about 0.47 in the previous analysis).
- sigma = k * RANDUNC: rescales every day equally, so the relative weights between days stay the same (the reason Gelman et al. give for rejecting the beta-binomial).

Implementation: the observation model is a config setting (`randunc_only` now; `randunc_plus_model_error` for the sensitivity run), so switching needs no code change.

### L2. Error distribution (Decided for now, stated limitation)

- Gaussian. Hollinger and Richardson (2005), the origin of RANDUNC, argue flux errors look more like a Laplace distribution. This is reported as a limitation.
- If a Laplace version is tried later, the prior on sigma must be reconsidered, because sigma means something different in the two distributions (Gelman et al. 2020, Section 7.1).

---

## Part 4. Findings so far

**F1. NEE is never missing in the daily file.** FLUXNET gap-fills every daily value, so the presence of NEE says nothing about measurement. Only QC does.

**F2. Fully gap-filled days.** 222 calibration days (4.3%) have QC = 0, meaning their NEE is entirely estimated.

**F3. QC and RANDUNC exclude different days** (calibration 1997 to 2010):

|  | RANDUNC missing or zero | RANDUNC OK |
|---|---|---|
| **QC below 0.75** | 116 | 335 |
| **QC at least 0.75** | 249 | **4,413 (assimilable)** |

- 335 days have a sigma but fail QC: RANDUNC alone would wrongly keep them.
- 249 days pass QC but have no sigma: days mostly filled with high confidence but with no truly measured half-hours.

**F4. Verification against the previous repository.** The new code reproduces its counts exactly: 4,413 of 5,113 calibration days (86.3%) and 810 of 1,461 evaluation days (55.4%).

**F5. Missing RANDUNC.** 1,138 of 6,940 days (16.4%) in the full file have no RANDUNC, including all of 2014. Within the calibration block, 365 days lack it, spread over several gaps (to be shown in Figure A), not one whole year.

**F6. Evaluation block.** 608 of 1,461 evaluation days (41.6%) are excluded only because RANDUNC is missing: all of 2014 and the first half of 2011.

---

## Part 5. Open items

**OPEN-1. Evaluation mask.** Evaluation needs sigma only for probability-based checks (for example, predictive interval coverage). Error measures like RMSE need only NEE and QC. Using QC alone keeps 1,418 instead of 810 evaluation days (adds 2014 and early 2011). Decide when designing the evaluation chapter.

### Carried over from the previous repository, not yet re-assessed

These will be reviewed when the work reaches them:

- Dirichlet reparameterisation of the allocation fractions.
- No EDC constraints (they break differentiability for NUTS).
- No VPD driver.
- NUTS as the only sampler.
- ACM photosynthesis implemented from Chuter et al. (2015), not Williams et al. (1997).
- Daily temperature extremes derived from the half-hourly file (the day/night proxy understated the range).
- LAI convention (hemisurface vs projected) kept as a two-way sensitivity run.
- `d_fall` inert as transcribed from Bloom and Williams (2015), equation A8.

---

## References

- Bacour, C. et al. (2023). Assimilation of multiple datasets results in large differences in regional- to global-scale NEE and GPP budgets simulated by a terrestrial biosphere model. Biogeosciences 20, 1089-1111.
- Douglas, et al. (2026). Regularisation of the 4DEnVar data assimilation method for calibration of land surface models. EGUsphere preprint.
- Famiglietti, C. A. et al. (2021). Optimal model complexity for terrestrial carbon cycle prediction. Biogeosciences 18, 2727.
- Fang, et al. Differentiable land model reveals global environmental controls on ecological parameters. arXiv:2411.09654.
- Fer, I. et al. (2018). Linking big models to big data: efficient ecosystem model calibration through Bayesian model emulation. Biogeosciences 15, 5801.
- Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808 (Section 10.5, pp. 53-54; Section 7.1).
- Hollinger, D. Y. and Richardson, A. D. (2005). Uncertainty in eddy covariance measurements and its application to physiological models. Tree Physiology 25.
- Kuppel, S. et al. (2012). Constraining a global ecosystem model with multi-site eddy-covariance data. Biogeosciences 9, 3757-3776.
- Pastorello, G. et al. (2020). The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Scientific Data 7, 225.
- Pu, et al. (2025). MCI GPP. Scientific Data.
- Richardson, A. D. et al. (2010). Estimating parameters of a forest ecosystem C model with measurements of stocks and fluxes as joint constraints. Oecologia 164, 25-40.
- Tramontana, G. et al. (2016). Predicting carbon dioxide and energy fluxes across global FLUXNET sites with regression algorithms. Biogeosciences 13, 4291.
- Williams, M. et al. (2009). Improving land surface models with FLUXNET data. Biogeosciences 6, 1341-1359.
- Yuan, Q. et al. (2025). Global carbon flux dataset generated by fusing remote sensing and multiple flux networks observation. Scientific Data.

Verify every citation against the original paper before using it in the thesis.