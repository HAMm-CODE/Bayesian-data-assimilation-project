# Decisions

Every methodological and technical choice in this project, with its reason.
Each entry is marked with a status:

- **Decided**: settled and implemented.
- **Proposed**: preferred option chosen, not yet final.
- **Open**: still to be decided; the options are recorded.

Last updated: 2026-10-12.

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
- Exception: Figure I highlights meteorological spring (March to May), set in the config as `eda.spring_months`.

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

### D16. Data and driver figures (Decided, 2026-10-08)

Script `scripts/02_data_eda.py`, separate from the sigma figures because it is about NEE and the drivers.

**Drivers shown: the ones DALEC2 reads (Figure H), plus shallow soil temperature in Figure I.** The model keeps its current drivers. A figure of every other driver in the file (VPD, precipitation, radiation components, wind, pressure, soil water and so on) was built and then removed by the author's choice, to keep the EDA focused. Adding a driver to the model would mean new equations, and VPD would reverse an inherited decision, so that stays future work.

Drivers DALEC2 reads: daily highest and lowest air temperature (from the half-hourly file, see D18) and their midpoint for respiration (D20), incoming shortwave radiation (`SW_IN_F`, converted from W m-2 to MJ m-2 d-1 by x 0.0864), atmospheric CO2 (`CO2_F_MDS`), day of year, and site latitude (61.85 N, fixed). Listed in the config under `drivers`; the full input list and the reason for each is D19.

The figures:

- **F. Data availability:** one row per year of the **whole record (1996 to 2014)**, one column per day, each day coloured as assimilable, QC below threshold, or QC passes but RANDUNC missing. Lines and labels mark the excluded year, the calibration block and the evaluation block. Whole record chosen so the figure shows why 1996 is excluded and why 2014 has no assimilable days.
- **G. NEE seasonal cycle:** one thin line per year, **calibration years and evaluation years in different colours**, with the calibration median on top. This replaces a separate calibration-versus-evaluation figure: it shows directly whether the evaluation years behave like the calibration years. Uses **QC-passing days** (sigma is not needed to describe NEE's shape; requiring it would remove 2014 entirely). Each line is a 15-day centred moving average (`eda.smoothing_days`).
- **Spring onset (in G):** first day of **10 days in a row with net uptake (NEE < 0)** on QC-passing days (`eda.onset_run_days`), computed per year. A day that fails QC breaks the run, so a data gap can delay the onset but never create one. Same definition as the previous repository. G shows the median onset of the calibration years and their range.
- **H. Model drivers:** air temperature (daily highest, mean and lowest; medians per day of year), incoming sunlight (median and middle half, plus potential sunlight with no atmosphere, which shows cloudiness), and CO2 as a time series (its upward trend is larger than its seasonal swing).
- **I. NEE against drivers:** NEE against incoming sunlight, air temperature and shallowest soil temperature, **spring days highlighted**. Spring is where light and temperature separate (light is already strong while air and soil are cold), so it tests which one limits uptake, and it is where DALEC2's onset is wrong (Stettz et al. 2022; Suni et al. 2003; Makela et al. 2004). Uses QC-passing days of the calibration block.
- **K. Day-to-day dependence:** autocorrelation of NEE **anomalies** (NEE minus the typical NEE for that day of the year) for gaps of 1 to 30 days, on assimilable calibration days. Removing the seasonal cycle first is necessary; otherwise "summer follows summer" looks like day-to-day dependence. Tests the likelihood's assumption that days are independent. The grey band shows the range expected if days were independent.

### D17. Shared plotting module (Decided, 2026-10-08)

- Figure style, colours, the 365-day calendar (29 February dropped), the per-day climatology (median and middle half) and the wrap-around smoothing live in `src/dalec/plotting.py`, used by both EDA scripts. One definition, so the two scripts cannot drift apart.

### D18. Daily highest and lowest temperature from the half-hourly file (Decided, 2026-10-12; closes OPEN-4)

**Why.** DALEC2's photosynthesis part (ACM) uses the daily temperature range (highest minus lowest). The daily FLUXNET file has no daily highest or lowest temperature. The daytime and night-time means (`TA_F_DAY`, `TA_F_NIGHT`) are averages over many hours, so they sit closer together than the true extremes; the previous repository found they understate the range by about 4.8 C. So the extremes are taken from the half-hourly file, as the previous repository did.

**How** (`daily_temperature_extremes` and `attach_temperature_extremes` in `src/dalec/data_io.py`; settings under `drivers.temperature_extremes` in the config):

- Read only two columns of the half-hourly file: `TIMESTAMP_START` and `TA_F` (gap-filled air temperature).
- Each half-hour belongs to the calendar day of its start time. Same rule as FLUXNET's own daily file, so the days line up.
- Daily highest = largest of the 48 half-hours; daily lowest = smallest. Stored as `t_max` and `t_min`.
- A day needs all 48 half-hours. With fewer, the true highest or lowest might be the missing one, so the day gets NaN instead of a partial answer. The count of valid half-hours is kept (`t_half_hours_valid`) and Figure H prints how many days are short.
- The run stops (fail loudly) if: a half-hour appears twice; a day has more than 48 half-hours; a date of the daily file is missing from the half-hourly file; on full days the mean of the half-hours differs from the daily `TA_F` by more than 0.05 C; or the daily mean lies outside the lowest-to-highest range. The mean check proves both files are the same version and the dates are not shifted.
- `TA_F_DAY` and `TA_F_NIGHT` are no longer used.

**Tested** on made-up data: correct extremes, one incomplete day set to NaN, and every check above raises an error when it should. Still to confirm on the real files.

### D19. Model inputs: only what DALEC2 reads (Decided, 2026-10-12)

**Rule.** The model receives only the inputs its equations use. Every extra column would either do nothing or need new equations and new parameters, which makes the parameters harder to pin down (identifiability).

| Input | Source | Where DALEC2 uses it |
|---|---|---|
| Respiration temperature | (highest + lowest) / 2, from the two rows below (D20) | Speed of decomposition and respiration (`exp(Theta * T)`) |
| Daily highest air temperature | half-hourly `TA_F` (D18) | Photosynthesis (ACM `exp(a8 * Tmax)` term); frost check |
| Daily lowest air temperature | half-hourly `TA_F` (D18) | Temperature range in canopy conductance (ACM); frost check |
| Incoming sunlight | `SW_IN_F`, daily file, x 0.0864 to MJ m-2 d-1 | Light supply to photosynthesis (ACM) |
| Atmospheric CO2 | `CO2_F_MDS`, daily file | Photosynthesis (ACM) |
| Day of year | from the date, not a column | Leaf-on and leaf-off timing (phenology); day length |
| Latitude 61.85 N | fixed site constant | Day length |

**Not inputs, on purpose:**
- `VPD_F` (air dryness) and `P_F` (precipitation). DALEC2 (Bloom and Williams 2015) has no water cycle; its ACM uses a fixed leaf water potential instead of a weather input. Adding them would mean new equations, so they stay future work. The old repository made the same choice ("No VPD"). Note: `research_plan_v2.md` still lists `VPD_F` and should be corrected.
- `TA_F_DAY` and `TA_F_NIGHT`: replaced by true extremes (D18).
- Daily `TA_F`: read only to check the half-hourly file (D18), never passed to the model (D20).
- `SW_IN_POT` and soil temperature (`TS_F_MDS_1`): shown in Figures H and I for context only. They sit under `eda.figure_only_columns` in the config, so they cannot be mistaken for model inputs.

**Checked against the literature (2026-10-12).** Bloom and Williams (2015, Appendix A) drive DALEC2 with "day of year, atmospheric CO2, minimum and maximum temperature, and global radiation", and state that "hydrological processes are not explicitly portrayed in the model". REFLEX (Fox et al. 2009, Table 1) uses the same four weather drivers, with the water terms of ACM "given a fixed value". VPD, wind speed and precipitation appear only in the later water-coupled ACM-GPP-ET (Smallman and Williams 2019), which also computes evaporation and is outside these research questions.

### D20. Respiration temperature = midpoint of daily highest and lowest (Decided, 2026-10-12)

**Choice.** The temperature in DALEC2's decomposition and respiration term is `t_mid = (t_max + t_min) / 2`, not the 24-hour mean `TA_F`.

**Why.** Bloom and Williams (2015) list only minimum and maximum temperature among DALEC2's drivers (Appendix A), so any daily temperature the model uses must come from those two. Chuter et al. (2015), whose written-out DALEC equations are the ACM source this project uses, state it directly: "Tm(t) is the mean daily air temperature, which is the mean of the maximum and minimum daily temperatures."

**What the literature does NOT say (checked 2026-10-12, correcting an earlier note).** The convention is not universal or explicit:
- Bloom and Williams (2015) use `T_t` in equations A5, A6 and A10 but never define it.
- Pinnington et al. (2016), DALEC2 at the Alice Holt flux site, say only "the daily mean temperature", without saying how it is computed.
- Delahaies et al. (2017) do not define the daily temperature either (an earlier version of this entry wrongly credited them with the Chuter quote).
- The previous repository used the 24-hour mean `TA_F` for respiration. D20 deliberately changes that.

**Option not taken.** `TA_F`, the true average of 48 half-hours, is slightly more accurate physically and is a defensible reading of "daily mean temperature". It was not chosen because the only explicit definition found (Chuter et al. 2015) is the midpoint, and the DALEC2 driver list contains only minimum and maximum.

**Size of the difference.** Figure H prints the average and spread of `t_mid - TA_F` on the real data, so the thesis can state how much this choice matters.

**Daily `TA_F` is still read**, but only for the D18 safety check (`drivers.temperature_extremes.daily_check_column`).

**Possible sensitivity check.** Because the two readings are both defensible, a short run with `TA_F` instead of `t_mid` would show whether the choice moves the posterior. Worth doing only if Figure H shows a sizeable difference.

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

**F7. Daily temperature extremes are not in the daily file (resolved by D18).** DALEC2's photosynthesis (ACM) uses the daily temperature range, but the FLUXNET daily file has no daily maximum or minimum. `TA_F_DAY` and `TA_F_NIGHT` are only stand-ins; the previous repository found they understate the range by about 4.8 C and derived true extremes from the half-hourly file. This is the one model input the daily file cannot supply directly.

**F9. Temperature-related cautions from earlier DALEC users (literature, 2026-10-12).**
- Bloom and Williams (2015) chose test sites "with a ≤ 3 months of recorded below-freezing soil temperatures" and little water stress, because "these criteria reflect the current capabilities of DALEC2". FI-Hyy, a boreal site, has a longer frozen season (the previous repository's notes say about four to six months of soil frost; to be confirmed from the soil temperature data). Using DALEC2 here is a deliberate scope choice and a stated limitation.
- Stettz et al. (2022), DALEC2 at a subalpine evergreen forest (US-NR1), report that baseline DALEC2 has no cold-temperature control on spring photosynthesis: GPP is "limited primarily by incoming shortwave radiation". They added a scalar on daily minimum temperature (GPP gradually reduced below about 6.0 C, fully stopped below about -7.1 C), which improved spring GPP (r = 0.88 vs 0.47). This matches the early spring onset found in the previous repository and supports showing spring in Figure I.

**F8. Driver completeness.** In the full record, the model drivers are complete except CO2 (0.3% empty, expected to be the 19 days at the start of 1996). Soil temperature and soil water are 1 to 5% empty; net radiation 21%, diffuse shortwave 46%, soil heat flux 46%, diffuse photosynthetic light 77%.

**F6. Evaluation block.** 608 of 1,461 evaluation days (41.6%) are excluded only because RANDUNC is missing: all of 2014 and the first half of 2011.

---

## Part 5. Open items

**OPEN-1. Evaluation mask.** Evaluation needs sigma only for probability-based checks (for example, predictive interval coverage). Error measures like RMSE need only NEE and QC. Using QC alone keeps 1,418 instead of 810 evaluation days (adds 2014 and early 2011). Decide when designing the evaluation chapter.

**OPEN-3. Spring onset number.** The previous EDA found a median observed onset of day 96 (range 58 to 117) with the 10-day rule, while project notes mention an observed onset of day 110 in the model comparison. These probably come from different definitions or data. Pick one definition (D16 uses the 10-day rule) and use it consistently in the thesis.

**OPEN-4. Daily temperature extremes (F7).** Closed 2026-10-12: derived from the half-hourly file (D18).

### Carried over from the previous repository, not yet re-assessed

These will be reviewed when the work reaches them:

- Dirichlet reparameterisation of the allocation fractions.
- No EDC constraints (they break differentiability for NUTS).
- No VPD driver.
- NUTS as the only sampler.
- ACM photosynthesis implemented from Chuter et al. (2015), not Williams et al. (1997).
- LAI convention (hemisurface vs projected) kept as a two-way sensitivity run.
- `d_fall` inert as transcribed from Bloom and Williams (2015), equation A8.

---

## References

- Bacour, C. et al. (2023). Assimilation of multiple datasets results in large differences in regional- to global-scale NEE and GPP budgets simulated by a terrestrial biosphere model. Biogeosciences 20, 1089-1111.
- Bloom, A. A. and Williams, M. (2015). Constraining ecosystem carbon dynamics in a data-limited world: integrating ecological "common sense" in a model-data fusion framework. Biogeosciences 12, 1299-1315.
- Chuter, A. M., Aston, P. J., Skeldon, A. C. and Roulstone, I. (2015). A dynamical systems analysis of the data assimilation linked ecosystem carbon (DALEC) models. Chaos 25(3), 036401. (Journal details from the previous repository; the quote in D20 was read from the manuscript at personalpages.surrey.ac.uk/a.skeldon/DALEC_paper.pdf, dated 8 December 2014. Check against the published article.)
- Delahaies, S., Roulstone, I. and Nichols, N. (2017). Constraining DALEC v2 using multiple data streams and ecological constraints: analysis and application. Geoscientific Model Development 10, 2635-2650.
- Douglas, et al. (2026). Regularisation of the 4DEnVar data assimilation method for calibration of land surface models. EGUsphere preprint.
- Famiglietti, C. A. et al. (2021). Optimal model complexity for terrestrial carbon cycle prediction. Biogeosciences 18, 2727.
- Fang, et al. Differentiable land model reveals global environmental controls on ecological parameters. arXiv:2411.09654.
- Fer, I. et al. (2018). Linking big models to big data: efficient ecosystem model calibration through Bayesian model emulation. Biogeosciences 15, 5801.
- Fox, A. et al. (2009). The REFLEX project: comparing different algorithms and implementations for the inversion of a terrestrial ecosystem model against eddy covariance data. Agricultural and Forest Meteorology 149, 1597-1615.
- Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808 (Section 10.5, pp. 53-54; Section 7.1).
- Hollinger, D. Y. and Richardson, A. D. (2005). Uncertainty in eddy covariance measurements and its application to physiological models. Tree Physiology 25.
- Kuppel, S. et al. (2012). Constraining a global ecosystem model with multi-site eddy-covariance data. Biogeosciences 9, 3757-3776.
- Pastorello, G. et al. (2020). The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Scientific Data 7, 225.
- Pu, et al. (2025). MCI GPP. Scientific Data.
- Pinnington, E. M. et al. (2016). Investigating the role of prior and observation error correlations in improving a model forecast of forest carbon balance using Four-dimensional Variational data assimilation. Agricultural and Forest Meteorology 228-229, 299-314.
- Richardson, A. D. et al. (2010). Estimating parameters of a forest ecosystem C model with measurements of stocks and fluxes as joint constraints. Oecologia 164, 25-40.
- Smallman, T. L. and Williams, M. (2019). Description and validation of an intermediate complexity model for ecosystem photosynthesis and evapotranspiration: ACM-GPP-ETv1. Geoscientific Model Development 12, 2227-2253.
- Stettz, S. G. et al. (2022). Resolving temperature limitation on spring productivity in an evergreen conifer forest using a model-data fusion framework. Biogeosciences 19, 541-558.
- Tramontana, G. et al. (2016). Predicting carbon dioxide and energy fluxes across global FLUXNET sites with regression algorithms. Biogeosciences 13, 4291.
- Williams, M. et al. (2009). Improving land surface models with FLUXNET data. Biogeosciences 6, 1341-1359.
- Yuan, Q. et al. (2025). Global carbon flux dataset generated by fusing remote sensing and multiple flux networks observation. Scientific Data.

Verify every citation against the original paper before using it in the thesis.