# Decisions

Every methodological and technical choice in this project, with its reason.
Each entry is marked with a status:

- **Decided**: settled and implemented.
- **Proposed**: preferred option chosen, not yet final.
- **Open**: still to be decided; the options are recorded.
- **Replaced**: a later entry changed it. The old text stays for the record, with a note naming the entry that replaced it.

Last updated: 2026-10-13 (full check of every entry against the current code, config and results).

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
- **Added for the model (2026-10-13):** PyMC and PyTensor (the gradient model, D46, and the sampler), pytest (runs the tests in `tests/`), and netcdf4 (saves the posterior draws as `.nc` files; ArviZ needs it to write them). Tested on the cloud test machine with PyMC 6.3.2 and PyTensor 3.3.3. `pymc`, `pytest` and `netcdf4` go into `environment.yml` from conda-forge like everything else; PyTensor is not listed on its own, because PyMC installs the PyTensor version it needs, and naming it separately could ask for a version that does not fit; the versions on the author's machine are the ones frozen in the lock file.
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

### D14. Config-driven loaders (Decided; updated 2026-10-13)

- Scripts load data only through the loaders in `src/dalec/data_io.py`, which read every setting from the config. There are three, one per need:
  - `load_screened_block()`: NEE, sigma, QC and the masks for one block (script 01).
  - `load_record_with_drivers()`: the whole record with the masks, the model drivers and the daily highest, lowest and midpoint temperature (script 02).
  - `load_model_block()`: the same for one block, which is what the model runs on (script 03, and the fitting scripts to come). Added in D41.
- (Earlier text said "a single entry point, `load_screened_block()`". That stopped being true once the model needed the drivers.)
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

Drivers DALEC2 reads: daily highest and lowest air temperature (from the half-hourly file, see D18) and their midpoint for respiration (D20), incoming shortwave radiation (`SW_IN_F`, converted from W m-2 to MJ m-2 d-1 by x 0.0864), atmospheric CO2 (`CO2_F_MDS`), the date (as day of year for day length, and as a running day count for leaf timing, D30), and site latitude (61.8474 N, fixed). Listed in the config under `drivers`; the full input list and the reason for each is D19.

The figures:

- **F. Data availability:** one row per year of the **whole record (1996 to 2014)**, one column per day, each day coloured as assimilable, QC below threshold, or QC passes but RANDUNC missing. Lines and labels mark the excluded year, the calibration block and the evaluation block. Whole record chosen so the figure shows why 1996 is excluded and why 2014 has no assimilable days.
- **G. NEE seasonal cycle:** one thin line per year, **calibration years and evaluation years in different colours**, with the calibration median on top. This replaces a separate calibration-versus-evaluation figure: it shows directly whether the evaluation years behave like the calibration years. Uses **QC-passing days** (sigma is not needed to describe NEE's shape; requiring it would remove 2014 entirely). Each line is a 15-day centred moving average (`eda.smoothing_days`).
- **Spring onset (in G):** first day of **10 days in a row with net uptake (NEE < 0)** on QC-passing days (`eda.onset_run_days`), computed per year. A day that fails QC breaks the run, so a data gap can delay the onset but never create one. Same definition as the previous repository. G shows the median onset of the calibration years and their range.
- **H. Model drivers:** air temperature (daily highest, midpoint and lowest, the three values the model uses; medians per day of year), incoming sunlight (median and middle half, plus potential sunlight with no atmosphere, which shows cloudiness), and CO2 as a time series (its upward trend is larger than its seasonal swing).
- **I. NEE against drivers:** NEE against incoming sunlight, air temperature (the midpoint, D20) and shallowest soil temperature, **spring days highlighted**. Spring is where light and temperature separate (light is already strong while air and soil are cold), so it tests which one limits uptake, and it is where DALEC2's onset is wrong (Stettz et al. 2022; Suni et al. 2003; Makela et al. 2004). Uses QC-passing days of the calibration block.
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

**Tested** on made-up data: correct extremes, one incomplete day set to NaN, and every check above raises an error when it should. **Confirmed on the real files (2026-10-12, F10):** every day has all 48 half-hours and every check passed.

### D19. Model inputs: only what DALEC2 reads (Decided, 2026-10-12)

**Rule.** The model receives only the inputs its equations use. Every extra column would either do nothing or need new equations and new parameters, which makes the parameters harder to pin down (identifiability).

| Input | Source | Where DALEC2 uses it |
|---|---|---|
| Respiration temperature | (highest + lowest) / 2, from the two rows below (D20) | Speed of decomposition and respiration (`exp(Theta * T)`) |
| Daily highest air temperature | half-hourly `TA_F` (D18) | Photosynthesis (ACM `exp(a8 * Tmax)` term) |
| Daily lowest air temperature | half-hourly `TA_F` (D18) | Temperature range in canopy conductance (ACM) |
| Incoming sunlight | `SW_IN_F`, daily file, x 0.0864 to MJ m-2 d-1 | Light supply to photosynthesis (ACM) |
| Atmospheric CO2 | `CO2_F_MDS`, daily file | Photosynthesis (ACM) |
| Day of year | from the date, not a column | Day length only (D23) |
| Day count `t` | from the date: days since the start of the run, day 1 = 1 January 1997 | Leaf-on and leaf-off timing (phenology, D30) |
| Latitude 61.8474 N | fixed site constant (`site.latitude_deg`) | Day length |

The midpoint temperature is also what the optional frost switch tests (D24); the switch is off in the main model.

**Not inputs, on purpose:**
- `VPD_F` (air dryness) and `P_F` (precipitation). DALEC2 (Bloom and Williams 2015) has no water cycle; its ACM uses a fixed leaf water potential instead of a weather input. Adding them would mean new equations, so they stay future work. The old repository made the same choice ("No VPD"); this re-assesses and keeps it. Note: `research_plan_v2.md` still lists `VPD_F` and should be corrected (OPEN-6).
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

**Size of the difference (real data, calibration block, 2026-10-12).** `t_mid - TA_F` averages -0.10 C; 90% of days lie between -1.05 and +0.78 C. Thesis sentence: the midpoint differs from the 24-hour mean by 0.1 C on average.

**Confirmed by the official code (2026-10-13).** The DALEC2 source code (CARDAMOM, `C/projects/CARDAMOM_MODELS/DALEC/DALEC_400/DALEC_400.c`, "Code used by Bloom et al., 2016. See also Bloom & Williams 2015") computes the decomposition multiplier as `exp(pars[9]*0.5*(maxtemp + mintemp))`, i.e. the midpoint. D20 matches the model as its authors run it, not only one paper's wording.

**Daily `TA_F` is still read**, but only for the D18 safety check (`drivers.temperature_extremes.daily_check_column`).

**Possible sensitivity check.** Because the two readings are both defensible, a short run with `TA_F` instead of `t_mid` would show whether the choice moves the posterior. Worth doing only if Figure H shows a sizeable difference. **Result: not needed.** An average difference of 0.1 C is small next to the 6.3 C daily range and the seasonal swing, so the run is dropped.

---

## Part 2b. Model structure (decided 2026-10-13; details in MODEL.md)

**Reference used.** Besides the papers, the official DALEC2 code (CARDAMOM, `DALEC_400`) was read: `DALEC_400.c`, `PARS_INFO_400.c`, `DALEC_ALL/ACM.c`, `DALEC_ALL/offset.c` (github.com/CARDAMOM-framework/CARDAMOM). Its header says it is the code used by Bloom et al. (2016), with Bloom and Williams (2015) as the model description. Where paper and code disagree, the code shows what the authors actually run. Caveat: its prior ranges are the 2016 ones and differ from the 2015 Table 1 (relevant to priors, not structure).

### D21. ACM in the Chuter et al. (2015) form (M1) (Decided)
Day length from latitude and solar declination, conductance `|psi_d|^a10 / (0.5 Tr + a6 R_tot)`, temperature on `Tmax`. **Why:** `ACM.c` in the official code computes exactly this. The Williams et al. (1997) form has no latitude in its day length and would overstate winter GPP at 62 N. This was an assumption in the previous repository; it is now confirmed.

### D22. ACM coefficients are DALEC2's fixed constants, at full precision (M2) (Decided)
`a2 = 0.0156935, a3 = 4.22273, a4 = 208.868, a5 = 0.0453194, a6 = 0.37836, a7 = 7.19298, a8 = 0.011136, a9 = 2.1001, a10 = 0.789798, psi_d = -2.0 MPa, R_tot = 1.0` (`DALEC_400.c`, `constants[]` and `gpppars[]`). **Why:** these are not a site choice; DALEC2 fixes them. They equal the "Loobos" set of Chuter et al. (2015), which the previous repository used in rounded form (for example 0.0156). They were fitted outside the boreal zone; stated limitation.

### D23. Solar day counted from 21 December (`doy + 10`) (M3) (Decided)
**Why:** the declination in `ACM.c` uses `(yearday + 10)`, matching Chuter et al. Section II.A. Without it the seasonal cycle of day length is 10 days late.

### D24. No frost cut-off in the main model; kept as a switch (M4) (Decided)
The previous repository stopped GPP when `(Tmax + Tmin)/2 < -2 C`. **Why not by default:** no such rule exists in the official code, and the cited precedent (Williams et al. 1997, Oregon conifers) has not been verified. Without it the thesis can state that DALEC2 is calibrated as distributed. (Note: the code's `gpppars[8] = -2.0` is the leaf water potential `psi_d`, not a frost threshold.) **Kept as a switch** because winter days carry most of the likelihood weight and DALEC2 lacks cold limitation on photosynthesis (F9, Stettz et al. 2022); it is tested once in the prior predictive check and once as a sensitivity run. **Status:** the prior predictive part is done (reported alongside in F12 to F15, never used to choose, D40); the sensitivity fit is still to do (planned runs list, Part 5). Config: `model.frost_cutoff` (`enabled: false`, threshold -2.0 C on the midpoint temperature).

### D25. A8 corrected: leaf fall anchored on `d_fall` (M5) (Decided)
`Phi_fall` uses `t - d_fall + psi_f`, not the printed `t - cr_fall + psi_f`. **Why:** the official code anchors leaf fall on its "Fday" parameter (`pars[14]`), the leaf-fall day, mirroring how leaf onset is anchored on "Bday". The printed A8 is a typo. This also removes the previous repository's inert `d_fall`.

### D26. All five phenology parameters are fitted, with site-based priors (M6) (Decided)
`d_onset`, `cr_onset`, `d_fall`, `cr_fall`, `c_lf` are sampled. **Why:** RQ1 and RQ3 ask what NEE can constrain; fixing parameters would hide the answer. Their priors were set in step 3 from Scots pine knowledge at Hyytiala (D37, with `c_lf` narrowed in D43), not from the 1 to 365 day published range.

### D27. Starting pools are sampled, with priors centred on site measurements (M7) (Decided)
All six `C(0)` are parameters with priors from Hyytiala stock measurements (Ilvesniemi et al. 2009; set in step 3). **Why:** the simplest first model (Gelman et al. 2020, Section 2.1). **Fallback:** the previous repository's reparameterisation (sample soil respiration and its split, derive the litter and soil pools) if the sampler struggles with the strong link between a turnover rate and its pool.

### D28. Allocation through a Dirichlet split (M8) (Decided)
`f_auto` plus a Dirichlet split of the remaining `1 - f_auto` across labile, foliage, roots and wood. **Why:** the official code instead splits step by step (foliage share of NPP, then labile share of the rest, then roots, wood taking the remainder). Both forms reach exactly the same set of allocations, so the model is the same; only the way the prior is written differs. The Dirichlet makes site-based priors simple and never produces negative fractions. This also re-frames correction 1 in `MODEL.md`: the paper's footnote is a shorthand of the code's step-by-step split.

### D29. LAI convention moved to the priors (M9) (Decided; closed by D38)
It only changes the prior range of `lma`, so it is decided with that prior in step 3. **Closed by D38:** the `lma` prior is wide enough to cover both LAI conventions, so no separate choice is needed.

### D30. Phenology uses a continuous day counter (M10) (Decided)
The leaf onset and leaf fall pulses use `t`, the number of days since the start of the run, not the day of year. **Why:** the official code does this (`MET[:,0]`, "projday"). With a 365.25-day cycle, a day-of-year counter that restarts every 1 January causes a small jump each year. Day length (B5) still uses the day of year, as in the code.

### D31. Leaf-fall offset `psi` from the official code's polynomial (A9) (Decided)
`psi` is computed with the sixth-order polynomial in `log(L - 1)`, `L = 1 / c_lf`, from `DALEC_ALL/offset.c`. **How it was settled:** the polynomial matches the equation `2 sqrt(pi) psi = log(1 - c_lf) exp(-psi^2)` within 0.004 for leaf lifespans 1.001 to 8 years; it does not match the previous repository's form (`log(1 - c_lf)` on the other side), which gave an offset about 8 times too large at `c_lf = 0.25` (about 24 days too early or late for a 60-day fall period). The A7 constants 0.6245 and 6.9088 come out of the same polynomial at `L = 1.001`, which confirms the reading. **Why the polynomial, not an exact root:** `c_lf` is now sampled (D26), so `psi` must change smoothly with it for NUTS gradients; the polynomial is differentiable, the root solve is not without extra work, and it is what DALEC2's authors run. Limit: valid for `L` in 1.001 to 8 years only, which covers Scots pine needles (3 to 5 years). **Still to check:** the printed A9 in the PDF, to tell a paper typo from a copying error.

### D32. Parameter scales (Gelman et al. 2020, Section 2.3) (Decided; final list updated after step 3)
Each parameter's sampling scale follows from its range by one rule: **log** if upper / lower >= 10, **interval** (logit) if less or if the parameter is a day-of-year position, **simplex** for the allocation split. **Final scales after step 3** (the rule was re-applied to the chosen priors, as planned):

| Scale | Parameters | Prior form |
|---|---|---|
| log | `theta_roo`, `theta_woo`, `theta_lit`, `theta_som`, `Theta` (D36), `cr_onset`, `cr_fall` (D37), `lma`, `ceff`, all six starting pools | log-normal |
| interval (0 to 1) | `f_auto`, `h` (D35), `c_lf` | beta |
| plain, unbounded | `d_onset`, `d_fall` (D37: their window is far from day 0 and day 365, so no bounds are needed) | normal |
| simplex | the four allocation shares | Dirichlet |

`theta_min` has no scale of its own: it is computed from `theta_lit` and `h` (D35).

(The first version of this entry, before step 3, put `Theta`, `cr_fall`, `d_onset`, `d_fall` on the interval scale and `theta_min` on the log scale. D35 to D37 changed those.)

**Why:** NUTS moves best when every parameter is unbounded and changes on a similar scale; the published ranges span from a factor of 4 to a factor of 10,000. The papers give ranges, not scales, so the scale is a modelling choice; using one rule keeps it simple and explainable. In PyMC most of this is automatic: a log-normal prior is sampled on the log scale, a bounded prior on the interval scale.

### D33. Prior strategy: smooth, weakly informative priors (Gelman et al. 2020, Section 2; step 3) (Decided)
Every prior is a smooth distribution on the parameter's scale (D32), centred on a site measurement or literature value and wide enough that the published range stays possible but its extremes are unlikely. No flat priors with hard edges. **Why:** hard edges are walls for NUTS; the previous repository's unusable run had eight parameters pinned at their prior bounds. Smooth priors also state the belief honestly: values near the measurement are likely, extremes possible but rare. **Options not taken:** flat priors on narrowed ranges (the previous repository's approach), and a mix of flat and smooth. **Cost:** each prior needs a stated centre and width with a source; recorded parameter by parameter below.

### D34. Priors, group 1: allocation (Decided)
- `f_auto ~ Beta(25.5, 24.5)`: centre 0.51, sd 0.07, 95% 0.37 to 0.65; under 0.4% outside the published 0.3 to 0.7. **Source:** Hyytiala tree growth 505 of GPP about 1,030 g C m-2 yr-1 (Ilvesniemi et al. 2009, Fig. 6) gives 0.51, site range 0.37 to 0.60; Waring et al. (1998) NPP/GPP 0.47 +/- 0.04 (f_auto 0.53); DeLucia et al. (2007) NPP/GPP 0.23 to 0.83 across forests, lowest in boreal forests (f_auto above 0.5).
- `(f_lab, f_fol, f_roo, f_woo) = (1 - f_auto) x Dirichlet(13.06, 13.06, 15.26, 44.27)`: shares of growth 0.152, 0.152, 0.178, 0.517 from the measured needle litterfall 154 (split equally between labile and foliage), below-ground litter 90 and wood growth 261 g C m-2 yr-1 (Ilvesniemi et al. 2009, Fig. 6); total concentration 85.65 derived from the measured spread of litterfall and GPP (previous repository, section 11). As shares of GPP: medians 0.073, 0.073, 0.085, 0.251.
- **Why this informative:** the previous repository's flat split sent too much carbon to the needles and LAI ran away to 48 (site value about 3); Gelman et al. (2020, Section 7.3): "as the model becomes more complex, all of the priors need to become tighter".
- **Caveats:** the equal labile/foliage split is a judgement (the data constrain only their sum); tower GPP includes the understory but the growth flows are trees only.
- **Planned check (Gelman et al. 2020, Section 6.3):** after fitting, compare prior and posterior spread; if the posterior of `f_auto` is about as wide as its prior, refit once with sd 0.10 as a sensitivity run.

### D35. Priors, group 2: turnover rates (Decided)
All log-normal (log scale, D32). Medians and 95% ranges given as residence times (1 / rate); `theta_lit` and `theta_som` are rates at 0 C.

| Parameter | Median | 95% (residence) | Source |
|---|---|---|---|
| `theta_roo` | 1.78e-3 d-1 (1.5 yr) | 0.6 to 3.7 yr, sdlog 0.45 | Hyytiala root litter 90 / root carbon 238 g C m-2 = 2.6 yr (Ilvesniemi et al. 2009); European conifers about 0.9 to 1.1 yr-1 (Brunner et al. 2013) |
| `theta_woo` | 5.33e-5 d-1 (51 yr) | 24 to 110 yr, sdlog 0.39 | fast end: site steady-state value (wood 6,300 / growth 261 g C m-2 yr-1), an upper limit on the rate because the stand is still growing; slow end: published lower limit of the rate |
| `theta_lit` | 1.22e-3 d-1 (2.2 yr) | 1 to 5 yr, sdlog 0.41 | Yasso07 litter decomposition rates (Tuomi et al. 2009), as worked out by the previous repository |
| `theta_som` | 7.15e-5 d-1 (38 yr) | 24 to 61 yr, sdlog 0.24 | Hyytiala soil carbon 6,560 g C m-2 and decomposition flux 290 to 370 g C m-2 yr-1 (Ilvesniemi et al. 2009), derivation from the previous repository |

**`theta_min` is not given its own prior. It is computed from the share of decomposing litter that becomes soil carbon**, `h = theta_min / (theta_lit + theta_min)`, so `theta_min = theta_lit h / (1 - h)`, with `h ~ Beta(4, 16)` (mean 0.20, 95% 0.06 to 0.40). **Source:** RothC sends 1/(x+1) of decomposed material to soil pools, with `x = 1.67 (1.85 + 1.60 exp(-0.0786 clay%))`: 17% at 5% clay, 19% at 10%, 22% at 23% (Coleman et al. 2024, RothC model description). Hyytiala's soil is a low-clay podzol, so about 0.2. Implied `theta_min`: median 2.8e-4 d-1, 95% 6.3e-5 to 1.1e-3, inside the published range on more than 99.99% of draws. **Why (Gelman et al. 2020):** the first proposal used the whole published range (a factor of 1,000) because no direct source was found, which conflicts with "as the model becomes more complex, all of the priors need to become tighter" (Section 7.3); rewriting the parameter as an interpretable quantity (Section 2.3) lets real subject-matter knowledge set the prior (Sections 5.9, 7.3). The model is unchanged; only the way its prior is written.

**Fallback rule for any parameter with no site or literature source:** centre on the geometric middle of the published range and let the published range be the 95% interval. Used only if no better source is found.

**Way forward for groups 1 and 2 (Gelman et al. 2020):**
1. **Judge the priors together, not one by one** (Sections 2.4, 7.3: independent priors "may jointly be strongly informative"). The prior predictive check (`scripts/03_prior_predictive.py`) runs DALEC2 with draws from all priors at once. Its pass criteria were written down **before** it was run (D40), so the check cannot be tuned after seeing it. **Done:** four rounds, F12 to F15, closed by D45.
2. **Watch the soil carbon used twice:** in `theta_som` here and in the starting soil pool (D27). The check above shows whether the implied soil respiration is sensible; if not, the previous repository's reparameterisation (D27 fallback) is used. **Result so far:** the soil pool stays steady (median end / start 1.00 to 1.04 in every round). Heterotrophic respiration is higher than the site value (median about 425 to 490 against 290 to 370 g C m-2 yr-1) but inside the Tier-2 band. The fallback is not needed for now; it stays ready in case the sampler struggles with `theta_som` and `c_som_0` (D45, safeguard 2).
3. **After fitting, measure how much each prior mattered** (Section 6.3): compare prior and posterior spread for every parameter. Parameters whose posterior is about as wide as their prior were not informed by NEE; this is a direct result for RQ1 (expected for the slow `theta_woo` and `theta_som`). Only those get a sensitivity refit with a wider prior.

### D36. Priors, group 3: temperature exponent (Decided, 2026-10-13)
`Theta ~ LogNormal(median 0.0336, sdlog 0.35)`: 95% 0.017 to 0.067, i.e. Q10 = exp(10 Theta) of 1.18 to 1.95 (median Q10 1.40). 3.7% below and 0.7% above the published 0.018 to 0.08. **Source:** Mahecha et al. (2010): ecosystem respiration across 60 FLUXNET sites converges to Q10 = 1.4 +/- 0.1 against air temperature once seasonal confounding is removed; the previous repository's estimate from Hyytiala winter fluxes was Theta = 0.0366 (Q10 1.44), close to the median. **Why wider than Mahecha's +/- 0.1:** DALEC2's single exponent also absorbs effects the model lacks (soil frost, moisture), so a site-level value may differ from the global one. **Scale changed from interval to log** (D32 allows this after step 3): the prior is now positive-only and smooth without needing the published bounds, which it already respects on 96% of draws.

### D37. Priors, group 4: leaf timing (phenology) (Decided; `c_lf` row replaced by D43)
In DALEC2 terms at an evergreen pine site, the onset pulse is the **spring flush of new shoots and needles** and the fall pulse is the **autumn drop of the oldest needles**. Note: this is **not** the NEE spring onset of Figure G (median day 100), which is the recovery of photosynthesis in needles that already exist; DALEC2 cannot represent that recovery (F9), and its leaf onset should not be tuned to it.

| Parameter | Prior | 95% | Basis |
|---|---|---|---|
| `d_onset` | Normal(145, 15) | day 116 to 174 | Scots pine shoot elongation in southern Finland starts in about mid to late May; with `cr_onset` 40 the release peaks around day 170 (mid-June), during needle elongation. **To confirm** with Schiestl-Aalto et al. (2013, Annals of Botany) or the dissertation (Schiestl-Aalto 2017), which give timing in heat sums |
| `cr_onset` | LogNormal(median 40 d, sdlog 0.35) | 20 to 79 days | length of shoot and needle elongation, roughly May to July; inside the published 10 to 100 |
| `d_fall` | Normal(265, 15) | day 236 to 294 | day of maximum needle fall (D25: the paper says the offset places the maximum loss rate on `d_fall`); Scots pine needle fall in Finland mostly late August to October. **To confirm** from monthly litterfall data |
| `cr_fall` | LogNormal(median 60 d, sdlog 0.35) | 30 to 119 days | spread of needle fall over about two months; inside the published 20 to 150 |
| `c_lf` | ~~Beta(14.51, 39.24)~~ **replaced by D43: Beta(43.84, 122.84)** | was 0.16 to 0.40 (needle life 2.5 to 6.2 years); now 0.199 to 0.332 (3.0 to 5.0 years) | mean needle retention 3.4 to 4.2 years in southern Finland (Pensa and Jalkanen 1999). The first version was wider than this source; D43 narrowed it |

Day positions use a plain normal prior (smooth, no bounds needed because the window is far from day 0 and day 365); `cr_fall` moves from interval to log scale (D32 re-applied). Days are counted on the continuous counter (D30), with day 1 = 1 January 1997.

### D38. Priors, group 5: canopy (`lma`, `ceff`), and the LAI convention (closes D29) (Decided; `ceff` replaced by D42)
- **`lma ~ LogNormal(median 185 g C m-2, sdlog 0.25)`**, 95% 113 to 302. **Basis:** foliage carbon from needle litterfall / `c_lf` = 154 / 0.26, about 590 g C m-2 (Ilvesniemi et al. 2009), divided by site LAI. Kolari (2010) gives all-sided LAI 6.5 to 8; as hemisurface LAI (divide by 2) this is 3.2 to 4.0, giving `lma` about 165; as projected LAI (divide by 2.57) 2.5 to 3.1, giving about 210. **The LAI convention question (M9) is resolved by putting it inside the prior:** the median sits between the two conventions and the 95% range covers both. One run instead of two sensitivity runs, and the uncertainty is stated honestly (Gelman et al. 2020, Section 7.3: uncertainty belongs in the prior). Note: about a third of the prior lies above the official code's upper limit of 200 (from Kattge et al. 2011), but inside the 2015 paper's 400.
- **Replaced by D42 (round 2): the prior now used is `ceff ~ LogNormal(median 10, sdlog 0.354)`, 95% 5 to 20.** The first version, kept for the record: **`ceff ~ LogNormal(median 16.5, sdlog 0.41)`**, 95% 7.4 to 37. **Basis:** `ceff` replaces `p11 x N` (nitrogen-use efficiency x foliar nitrogen per leaf area). `p11 = 7.4` belongs to the same ACM fit as the fixed coefficients (D22, Chuter et al. 2015, Loobos); with needle N between about 1 and 5 g N m-2 leaf area this gives 7.4 to 37, which also contains the Loobos value 29.6 and most of the previous repository's REFLEX range 5 to 20. (The "to check: Hyytiala needle N" note no longer sets the prior after D42; it is kept only for the discussion of what `ceff` means, D42 "honest tension".)

### D39. Priors, group 6: starting pools on 1 January 1997 (implements D27) (Decided)
All log-normal, centred on Hyytiala measurements (Ilvesniemi et al. 2009, Fig. 6 and Table 4).

| Pool | Median | sdlog | 95% | Basis |
|---|---|---|---|---|
| `c_lab_0` | 77 | 0.50 | 29 to 205 | one year of labile allocation (`f_lab` x GPP, about 0.073 x 1,030), held over winter for the spring flush |
| `c_fol_0` | 590 | 0.30 | 328 to 1,062 | needle litterfall 154 / `c_lf` 0.26 |
| `c_roo_0` | 238 | 0.30 | 132 to 429 | measured fine-root biomass 476 g m-2 dry mass, halved for carbon |
| `c_woo_0` | 5,900 | 0.25 | 3,600 to 9,600 | tree carbon 6,800 minus foliage, labile and fine roots |
| `c_lit_0` | 450 | 0.50 | 170 to 1,200 | litter input about 260 g C m-2 yr-1 times a field residence of about 1.8 years (D35 litter rate, D35 share to soil) |
| `c_som_0` | 6,100 | 0.25 | 3,700 to 10,000 | soil carbon 6,560 minus the litter pool |

**Caveats:** the measurements are from the 2000s while the run starts in 1997 (the stand was younger, and a thinning followed in 2002); the widths are set to absorb this. The litter pool rests on a steady-state argument and is the least certain, hence its wider prior. The soil measurement also informs `theta_som` (D35); the prior predictive check tests whether the implied soil respiration is sensible.

### Summary of step 3 (closed by D45)
24 sampled quantities: `f_auto`, the 4-share Dirichlet (3 free, because the shares add up to 1), 4 turnover rates, `h`, `Theta`, 5 phenology, `lma`, `ceff`, 6 starting pools. That is 23 free dimensions for the sampler; `theta_min` is computed, not sampled. Every prior is smooth (D33) and has a stated source or reason.

**Priors used for fitting:** D34 to D39, with `ceff` from D42 and `c_lf` from D43, plus the soft foliage constraint D44. All in one list, `PRIORS` in `src/dalec/priors.py` (D41).

**Still marked to confirm:** leaf onset timing and needle fall timing (D37). Needle N no longer sets a prior (D42). These are listed in OPEN-6.

### D40. Prior predictive check: pass criteria, written before running (Gelman et al. 2020, Section 2.4) (Decided)
**Set-up.** 1,000 draws from all priors (D34 to D39) at once; each draw runs DALEC2 over the full calibration block (1997 to 2010, 5,113 days) with the real drivers; frost switch off (main) and once on (D24). No NEE observations are used: the check asks whether the priors describe a **believable boreal pine forest**, not whether they fit the tower data. The observed NEE climatology may be shown in a figure for context only; nothing is tuned to it.

**Tier 1: physically impossible.** A draw fails if any of these happens on any day:
- a pool becomes negative or not a number;
- LAI (`C_fol / lma`) above 12;
- daily GPP above 20 g C m-2 d-1;
- annual GPP above 3,000 g C m-2 yr-1.
**Target: at most 5% of draws fail Tier 1.** (The previous repository's published-range priors: 84% impossible.)

**Tier 2: plausible for this forest.** Broad bands, about half to double the independent site values (Ilvesniemi et al. 2009; Kolari 2010), averaged over 1997 to 2010:
| Quantity | Band | Site value |
|---|---|---|
| annual GPP | 500 to 2,200 g C m-2 yr-1 | 952 to 1,104 |
| annual heterotrophic respiration | 145 to 740 | 290 to 370 |
| annual NEE | -700 to +300 (sink or source allowed) | a sink |
| peak LAI | 1.5 to 8 | about 2.5 to 4 |
| foliage pool, end / start | 1/3 to 3 | (no runaway or collapse) |
| wood pool, end / start | at least 0.8 | growing stand |
| soil pool, end / start | 0.7 to 1.3 | soil changes slowly |
**Target: at least half of the draws pass every Tier-2 check;** the pass rate of each check is reported separately. The site GPP value is itself partly derived from tower data, so its band is deliberately wide.

**Why these targets:** Gelman et al. (2020) set no numbers; the principle is that the prior should rule out the impossible while still allowing more than the plausible. The 5% and 50% are judgements, stated here before the check is run so they cannot be adjusted after seeing the results.

**If the check fails:** find which priors drive the failures (compare the parameter values of failing and passing draws), change only those, and only with a stated source or reason; record every change here and re-run. The frost-switch run is reported alongside, not used to pick between them.

### D41. Implementation of the priors and the prior predictive check (Decided, 2026-10-13)
- **One list of priors** (`src/dalec/priors.py`, `PRIORS`), used by the prior predictive check now and by the PyMC model later, so they cannot disagree.
- **`c_lf` draws outside the polynomial's range** (needles living longer than 8 years, about 0.3% of the prior) are drawn again from the same prior. This trims a tiny tail; D31's polynomial is not valid there. (With D43's narrower prior this almost never happens, but the rule stays as a guard.)
- **Prior predictive settings** in the config (`prior_predictive`): 1,000 draws, fixed seed 20261013, run in batches of 250.
- **Diagnosis when a target is missed:** for every sampled quantity, the standardised difference between failing and passing draws on the prior's own scale (log, logit or plain); the largest ones point to the priors to revisit (D40, "if the check fails").
- **The loader for model runs** (`load_record_with_drivers`, `load_model_block`) moved into `data_io.py`, shared by scripts 02 and 03.
- **The D40 checks** (`TIER2_BANDS`, `summarise_batch`, `evaluate`) moved from script 03 into `src/dalec/plausibility.py` (2026-10-13), because the fake-data check (D48) uses them too. One definition, so the two scripts cannot drift apart (as D17).
- **Shared fit steps** (chain starts, the prior-against-posterior table, sampler health, the D49 computation checks, the interval figure) moved from script 04 into `src/dalec/posterior.py` (2026-10-15), so the fake-data check and the real fits are judged the same way. Config: `diagnostics` holds the D49 thresholds and `sampler.quick` the quick-run settings for both scripts.

### D42. Prior predictive round 2: `ceff` prior (Decided, 2026-10-13, after the re-run, F13)
Following F12 and D40 ("change only the priors linked to failure, with a stated source"), only `ceff` changes in this round.
- **New:** `ceff ~ LogNormal(median 10, sdlog 0.354)`, 95% 5 to 20. **Source:** the REFLEX DALEC model-data fusion experiment used a prior of 5 to 20 for the same ACM parameter ("nitrogen use efficiency parameter (a1) in ACM", Fox et al. 2009, Table 3; the previous repository cited it as Table 4). REFLEX includes an evergreen forest site.
- **Replaces** D38's `LogNormal(16.5, 0.41)` (95% 7.4 to 37), which was built from `p11 x N` with an assumed needle N of 1 to 5 g N m-2.
- **Honest tension:** a typical Scots pine needle N per area (about 1.2% N of dry mass on roughly 370 g dry mass m-2) would point to N near 4 and `ceff` near 30, i.e. higher, not lower. The prior predictive check says that, with DALEC2's fixed ACM coefficients (fitted at a temperate site, D22) and no cold limit (F9), such values make GPP implausibly high at Hyytiala. So `ceff` here acts as a site-level efficiency that also absorbs ACM's mismatch at a boreal site; it cannot be read as a physiological nitrogen value. This link between `ceff` and the fixed ACM coefficients is an identifiability point for RQ3.
- **Gelman et al. (2020, Section 2.4)** allow priors to be revised after a prior predictive check, because no observations enter it; the revision must rest on outside knowledge (here a published prior), and both rounds are reported.
- **Cold-limit model:** see D57, which returns `ceff` to D38 in that model only, because the reason given here (no cold limit) does not apply to it.
- **Options not taken:** (B) setting `ceff` so that the implied annual GPP matches the site's plausible range (prior elicitation on an observable, also Section 2.4; rejected for now because the site GPP comes partly from the same tower); (C) keeping D38 and accepting the failure.

### D43. Prior predictive round 3: `c_lf` prior narrowed to its source (Decided, 2026-10-13, after the re-run, F14)
- **New:** `c_lf ~ Beta(43.84, 122.84)`: mean 0.263, 95% 0.199 to 0.332, i.e. needle lifespan **3.0 to 5.0 years**.
- **Replaces** D37's `Beta(14.51, 39.24)` (lifespan 2.5 to 6.2 years).
- **Why:** `c_lf` is the only prior linked to both failure modes (F13). D37's prior was **wider than its own source**: Pensa and Jalkanen (1999) report mean needle retention of 3.4 to 4.2 years in southern Finland; 3 to 5 years covers that with a margin (the previous repository used the same 3 to 5 years). Narrowing a prior to what its source supports is outside knowledge, not tuning to the check (D40, Gelman et al. 2020, Section 2.4).
- **Not changed in this round:** `lma` (linked to Tier 1, but its width reflects the unresolved LAI convention, D38); kept for round 4 if needed.

### D44. Soft steady-state constraint on the foliage pool (Decided, 2026-10-13, after the re-run, F15)
**What.** The joint prior gets one extra smooth term: `log p += -0.5 (log Rs / log 2)^2`, where for the foliage pool `Rs = (total inputs / total outputs over the run) x (mean pool on each 1 January) / (pool at the start)`. A foliage pool in balance has `Rs = 1` and no penalty; a pool running away (`Rs` above 1) or collapsing (below 1) is made less likely, smoothly. Code: `steady_state_ratio` and `edc_log_weight` in `priors.py`.

**Source: the official DALEC2 code.** `DALEC_400/EDC2_400.c` implements "EDCs 7-12" exactly this way for every pool, with equilibrium factor `EQF = 2` by default: it rejects a draw when `|log Rs| > log(EQF)` and also stores the Gaussian version `-0.5 (log Rs / log EQF)^2` ("EDCPROB"). The growth-rate and decay-rate rules described as EDC 7 and EDC 8 in Bloom and Williams (2015) are commented out in that code version and replaced by this steady-state form. Where paper and code differ we follow the code (as in D21 to D31).

**Why soft, and why foliage only.**
- Hard accept/reject rules break NUTS (the carried-over reason for dropping EDCs). The Gaussian form is smooth, so NUTS can follow its gradient: it turns "no EDCs" into "**soft EDCs only**". This revises a carried-over decision and should be raised with the supervisor.
- Foliage only, because the prior predictive check points to the canopy (F12 to F14). The official code applies it to all six pools, but Bloom and Williams (2015) chose sites near steady state; Hyytiala is a growing, managed stand (thinned in 2002), so a steady-state rule on wood would penalise real growth.
- In Gelman et al. (2020, Section 7.3) terms, it is a **joint prior**: a statement about what the parameters do together, which the independent priors could not express.

**How it is tested.** The prior predictive script re-weights the same 1,000 draws by the constraint (importance re-weighting, Gelman et al. 2020, Section 6.3) and reports the D40 rates under the new joint prior, plus the effective number of draws. No new sampling is needed.

**Honest caveat.** The constraint and the Tier-2 "foliage end/start" check measure related things, so part of the improvement is expected by construction. The justification for adding it is the official code, not the check.

### D45. Step 3 closed: the joint prior is D34 to D44; the narrow Tier 1 miss is accepted (Decided, 2026-10-13)
**Decision.** Stop revising priors. The prior used for fitting is the independent priors D34 to D43 plus the soft foliage constraint D44.
**Why stop here.** Four rounds took impossible draws from 84% (previous repository) and 34.7% (round 1) to 6.4%, and Tier 2 now passes. The remaining 6.4% are mostly days with GPP above 20 g C m-2 d-1; such draws imply daily NEE far beyond anything observed (summer observations about -2 to -4), so the likelihood rules them out directly. Further narrowing has small returns and would push priors past what their sources support (Gelman et al. 2020, Section 7.3: priors should encode knowledge, not be tuned until a check passes).
**The miss is reported, not hidden:** D40's 5% target is kept as written; the thesis states 6.4% against 5%, with this reasoning.
**Safeguards for fitting:**
1. Start every NUTS chain from a prior draw that passes Tier 1, so no chain begins in a broken region.
2. Watch divergences and chains stuck at extreme GPP or LAI; if they appear, revisit `lma` (strongest remaining Tier 1 link after `ceff`) or extend D44 to the labile pool.
3. The frost switch stays a sensitivity run (D24), not a fix chosen by this check.

### D46. The PyTensor model for NUTS (`src/dalec/model.py`) (Decided, 2026-10-13)
- **Same model as `model_numpy.py`, verified:** daily NEE, GPP, end pools and the D44 ratio agree to a relative 1e-10 for three parameter sets and with the frost switch (`tests/test_model.py`).
- **Gradients verified:** for 7 parameters spread over the groups (`ceff`, `lma`, `temperature_exponent`, `theta_som`, `c_lf`, `d_fall`, `c_som_0`), PyTensor's gradient of total NEE equals a finite-difference slope of the numpy model within 0.01%.
- **Organisation for speed:** weather-only terms are computed once in numpy; leaf timing and the temperature multiplier for all days at once outside the loop; only photosynthesis and the pool updates run inside `pytensor.scan`, with parameters passed as explicit `non_sequences` (the previous repository measured a 12x slower gradient through closures).
- **Measured cost (cloud test machine, 14 years of synthetic weather):** about 40 ms for one value plus gradient of all 24 parameters, after a one-off compile of about 35 s. A NUTS draw needs roughly 50 to 250 such evaluations.
- **Leaf timing uses the D31 polynomial** written as a PyTensor expression, so `c_lf` stays differentiable (D26).
- **Tests:** ACM 14, numpy model 14, priors 10, PyTensor model 11 (the last group takes a while because each case compiles the model). With the 4 tests of the PyMC model (D47), 53 in total.
- The PyMC model that joins the priors, the D44 constraint and the likelihood is D47; the fake-data check is D48 and D49.

---

## Part 2c. Fitting (Gelman et al. 2020, Sections 3 and 4; decided 2026-10-13)

### D47. The PyMC model (`src/dalec/inference.py`) (Decided, 2026-10-13)
**What it is.** log posterior = log priors (D34 to D43) + the soft foliage term (D44) + the log likelihood of daily NEE on assimilable days (L1, D12). It adds no new equations: the priors come from the one list in `priors.py` (D41) and the forward model from `model.py` (D46).
- **Priors:** log-normal and beta priors are sampled on the log and logit scales (D32); the allocation split is a Dirichlet times `1 - f_auto` (D28); `theta_min` is computed from `theta_lit` and `h` (D35).
- **`c_lf`:** a beta prior cut to needle lives of 1.001 to 8 years, the same rule as the numpy draws (D41). PyMC still samples it on a smooth, unbounded scale, so the cut is not a wall for NUTS.
- **Soft constraint:** `pm.Potential` with `-0.5 (log Rs / log 2)^2` (D44).
- **Likelihood:** `NEE_t ~ Normal(model_t, RANDUNC_t)` on assimilable days; the sensitivity version `sqrt(RANDUNC_t^2 + sigma_m^2)` with `sigma_m ~ HalfNormal(1)` is the same code with one config switch (`likelihood.observation_model`, L1).
- **Checked (`tests/test_inference.py`):** at two different parameter sets, PyMC's log posterior equals the same sum computed independently with scipy and the numpy model, to 1e-9 relative, for both observation models.
- **A bug found and fixed while testing:** PyTensor stores a plain Python number that fits exactly in 32 bits (such as 145.0 or 25.5) as a 32-bit number, so some prior densities were computed to only about 7 digits. Every prior constant is now made 64-bit. The effect on sampling would have been tiny, but the test now shows the model is exact.

### D48. Fake-data check: design (Gelman et al. 2020, Section 4.1) (Decided, 2026-10-13)
**Question:** if we knew the right answer, would the fit find it? Only if it does can a real-data result be trusted.
- **The true forest:** the first of 200 prior draws (fixed seed) that passes every D40 check, Tier 1 and Tier 2. Believable, but not the prior centre, which would make the test too easy. Taking the first one that passes, rather than choosing, avoids picking a convenient truth.
- **Fake NEE:** DALEC2 run with the true values on the real 1997-2010 weather, plus random normal noise with the real RANDUNC as its size, on the real assimilable days only. So the fake data have exactly the real gaps, the real winter-heavy weights and the real amount of information.
- **Fitted with the main model** (`randunc_only`), because the fake noise is RANDUNC only.
- **Options not taken:** the prior median as truth (too easy); several truths, or simulation-based calibration (Gelman et al. 2020, Section 4.2: hundreds of fits, too expensive at hours per fit). One fit is the minimum Gelman et al. recommend; a second truth can be added if the first shows problems.
- **Quick run first** (`--quick`): 2 years, 2 short chains. It only shows that every step runs and gives a time estimate; nothing is judged ("fit fast, fail fast", Section 3.4).

### D49. Sampler settings and pass criteria, written before running (Decided, 2026-10-13)
**Settings** (config `sampler`): 4 chains, 1,000 warm-up and 1,000 kept draws each, `target_accept = 0.9` (smaller, safer steps than PyMC's default 0.8, because the prior check showed a narrow canopy region), fixed seed. Each chain starts from a different prior draw that passes Tier 1 (D45 safeguard 1), with PyMC's random "jitter" of the start switched off, because it can move a start by a factor of about 3 on the log scale and out of the plausible region.

**Pass criteria for the fake-data check** (config `fake_data`):
1. no divergences (a divergence = the sampler's step broke down, a sign it cannot explore part of the posterior);
2. R-hat below 1.01 for every quantity (R-hat compares the chains; 1 = they agree);
3. effective sample size at least 400, bulk and tail, for every quantity (how many independent draws the chains are worth);
4. the true value inside the 90% posterior interval for at least 20 of the 24 sampled quantities.

**Why these numbers.** 1 to 3 are the current standard thresholds (Vehtari et al. 2021). For 4: even a perfect fit misses about 1 in 10 at 90%, so about 2.4 misses are expected by chance; 20 of 24 allows for that. The quantities are not independent, so this is a judgement, stated before running.

**Also reported, not pass/fail (RQ1 preview):** for each quantity, posterior width / prior width on the prior's scale. Below 0.5 = strongly informed by NEE; 0.5 to 0.9 = partly; 0.9 or more = not informed (D35, way forward 3). The prior widths are from the independent priors, without the D44 term.

**If it fails:** look for the cause before changing anything (Gelman et al. 2020, Section 5): divergences point to geometry (then consider the D27 reparameterisation or a higher `target_accept`); a truth outside a narrow interval points to a bug or a non-identified pair; record every change here.

### D50. Where the fake-data check runs: Roihu, with the existing container (Decided, 2026-10-13)
- **Roihu, partition `small`**, account `project_2020170`: one node, 4 CPUs (one per chain), 8 GB memory, time limit 48 h (`small` allows 72 h; only the time used is billed). Time raised from 24 h after the first full job timed out; the run took 29.4 h (F16).
- **Memory, measured (2026-10-15):** Roihu bills a `small` job on whichever is larger, 0.75 BU per CPU-hour or 0.375 BU per GiB-hour of memory requested. The first full run asked for 32 GB but used 1.52 GB (`seff 2278889`), so it was billed on memory: 353 BU instead of about 88. Full jobs now ask for 8 GB (2 GB per CPU, where the two charges are equal; over 5 times the measured use) and quick jobs for 4 GB. CPU efficiency was 99%, so one CPU per chain is fully used and more CPUs would not speed up a chain. Job scripts: `jobs/fake_data_quick.sbatch`, `jobs/fake_data_full.sbatch`, shared set-up `jobs/roihu_setup.sh`.
- **Python environment:** the Tykky container built for the previous repository (`/projappl/project_2020170/dalec-env`; Python 3.13, PyTensor 3.2.4, C++ compiler inside). Our code is not installed into it; it is found through `PYTHONPATH=src`, and the job stops if `dalec` would be loaded from anywhere else. **Checked:** all 53 tests pass with PyTensor 3.2.4. The PyMC and ArviZ versions in the container are printed at the start of every job.
- **Why not a new container tonight:** building one takes time and the existing one already ran a model of the same size (3.9 h for 4 x 2,000 draws in the previous repository). A fresh container from `environment.yml` is the clean option before thesis results (D2, lock file).
- **Order:** the quick job first; the full job is submitted with `--dependency=afterok`, so it starts only if the quick job finished without an error.

### D51. The real-data fit (`scripts/05_fit_real_data.py`) (Decided, 2026-10-15)
**What.** The model (D47) and sampler settings (D49) that passed the fake-data check (F16), fitted to the measured `NEE_VUT_REF` on the 4,413 assimilable days of 1997-2010. Chain starts: physically possible prior draws (D45) from a different seed than the fake-data check. Job: `jobs/fit_full.sbatch` (48 h), after `jobs/fit_quick.sbatch`.

**Reading order, fixed before running.**
1. **Computation first:** D49 criteria 1 to 3 (no divergences, R-hat below 1.01, effective sample size at least 400). If any fails, nothing else is interpreted until the cause is found.
2. **Prior against posterior** for every parameter (Figure O, same table as F16 without a truth).
3. **The forest the posterior implies:** annual GPP, autotrophic, heterotrophic and ecosystem respiration, NEE, peak LAI and the pool changes, next to the independent site values (RQ3). Also the share of posterior draws outside each D40 band.
4. **Posterior predictive checks** (Gelman et al. 2020, Section 6.1): 400 posterior draws run forward, each with a fresh error of the size the likelihood assumes. Reported: the share of assimilable days whose observed NEE lies inside the 90% predictive interval (overall and by season), the average misfit by season, annual sums on assimilable days, the seasonal cycle (Figure P) and the monthly series (Figure Q). These are descriptive: nothing is tuned to them.

**Expectations written down before seeing the results** (so they can be confirmed or refuted honestly):
- **Spring misfit:** the model turns to net uptake too early in spring (F12, Stettz et al. 2022), so in March to May the observed NEE should sit above the model (observed minus model positive).
- **Coverage below 90%:** RANDUNC covers measurement error only, not model error (L1), and days are not independent (F11), so the predictive interval is likely too narrow, most clearly where the model is structurally wrong (spring).
- **Weakly identified pairs:** root pool and root turnover, soil pool and soil turnover, and the starting wood pool should be only partly or not informed (F16).

**Sensitivity runs** use the same script with options (`--frost`, `--observation-model randunc_plus_model_error`, `--mass-matrix dense`); each writes its own files, so several can run at the same time on Roihu.

### D52. Speed test: dense mass matrix on the fake data (Decided, 2026-10-15)
**Why.** The fake-data fit took 29.4 h because nearly every draw used the maximum 1,023 gradient evaluations at a step size of about 0.005 (F16). A "mass matrix" is the sampler's learned picture of the posterior's shape. The validated setting ("diagonal") learns one scale per parameter; a "dense" one also learns how parameters move together, which can straighten out the sloping ridges of the pool and rate pairs and allow longer steps. With 23 parameters a dense matrix is cheap. It changes only how the sampler moves, not the model or the answer.

**Test.** The same fake-data check (D48), same truth, seeds and starts, with `--mass-matrix dense` (PyMC `init="adapt_full"`), run on Roihu alongside the real-data fit.

**Adopt dense for later runs only if all three hold** (fixed before running): (1) it passes all four D49 criteria; (2) it takes at least 30% less time than 29.4 h, i.e. at most about 20.5 h; (3) for every quantity, its posterior median differs from the diagonal run's by less than 10% of the diagonal run's 90% interval width. Otherwise the diagonal setting stays. The main real-data fit (D51) uses the validated diagonal setting either way.

### D53. Second speed test: dense mass matrix with target_accept 0.95 (Decided, 2026-10-16; adopted for full runs after F23)
**Why.** F19: dense was 6 times faster and gave the same posterior, but had 13 divergences. A divergence marks a place where the step was too large for the curvature; raising `target_accept` from 0.9 to 0.95 makes the sampler choose smaller steps, which usually removes them at a modest cost in speed.
**Test.** The same fake-data check (D48) with `--mass-matrix dense --target-accept 0.95`.
**Adopt dense with target_accept 0.95 for later runs only if all three D52 criteria hold** (all four D49 criteria including no divergences; at most about 20.5 h; every posterior median within 10% of the diagonal run's 90% interval width). Otherwise the diagonal setting stays. Fixed before running.
**Code:** scripts 04 and 05 accept `--target-accept`; the value is added to the output names (e.g. `_dense_ta95`).
**Outcome (F23): all three criteria hold, so full runs now use a dense mass matrix with target_accept 0.95** (`config/default.yaml`). Quick runs keep diagonal with 0.9: 150 warm-up draws are too few to learn a dense mass matrix (25 x 25 numbers), and the quick runs F17 to F22 stay comparable. Since this change, output names carry a sampler tag (`_dense`, `_diagonal`, `_taNN`) only when the setting is given on the command line; results made before carry the older names.

### D54. A simpler model: holding the canopy parameters fixed (Decided as a test, 2026-10-16)
**Why.** With model error (L1), the chains agree, but the fit still builds GPP from an implausible canopy: thin, partly deciduous, very efficient leaves (F18). Gelman et al. (2020, Section 5.2) suggest simplifying a poorly behaving model, "stripping it down until you get something that works". The aim is to take away the canopy's freedom to imitate a winter shutdown, so the remaining misfit shows up openly (in the predictive checks and in sigma_m) instead of hiding in the parameters.

**How.** `--fix name,...` in script 05 holds the named parameters at their prior medians (`prior_median` in `priors.py`), which are site-based (D37, D38, D42, D43). They are stored as constants, so every later step reads them like sampled ones; Figure O marks them "fixed". Tested: the log density equals the full one minus the fixed parameters' prior terms (`tests/test_inference.py`).

**The ladder of simpler models, from least to most fixed:**

| Version | Held fixed (value) | What it removes |
|---|---|---|
| S1 | `lma` (185 g C m-2), `c_lf` (0.262, needles about 3.8 years) | the canopy's freedom to become thin and short-lived |
| S2 | S1 + `ceff` (10, REFLEX median) | also the freedom to make leaves very efficient; GPP size then follows the weather and the fixed canopy |
| S3 | S1 + leaf calendar `d_onset`, `cr_onset`, `d_fall`, `cr_fall` | also the timing of needle flush and fall |
| S4 | S1 + the starting pools NEE cannot inform: `c_woo_0` (5,900), `c_roo_0` (238), `c_som_0` (6,100), all from site measurements (D39) | the pool/rate trade-offs of F16: with the pool fixed at its measured size, NEE can inform the turnover rate |

S1, S2 and S4 are tested first (quick runs, real NEE 1997-1998, model error on). S4 was added on 2026-10-16 at the author's request. Note: `research_plan_v2.md` already describes the model as "a reduced DALEC model", so a simpler version is in line with the original plan.

**Other simplifications considered:**
- *A structurally simpler DALEC:* the evergreen DALEC of Williams et al. (2005), with no labile pool and no leaf-timing pulses, only a constant needle loss rate. It cannot build a seasonal canopy at all. (An evergreen DALEC comparison was earlier dropped as an RQ3 plan; here it would be a simplification of the model, not a comparison.) Bigger change: new model code and new checks. **Not chosen** (author, 2026-10-16): DALEC2 stays as the model structure.

**The expansion alternative, for comparison (not chosen yet):** DALEC2cold (Stettz et al. 2022) multiplies GPP by a cold scalar on minimum temperature: 0 below T0, 1 above Tg, a straight line in between, with T0 and Tg calibrated (at Niwot Ridge T0 = -7.1 +/- 1.1 C, Tg = 6.0 +/- 2.6 C). Two extra parameters. Stettz et al. had no LAI constraint and note that their site's LAI is roughly constant through the season; they did not test whether baseline DALEC2 compensated through the canopy.

**How each version is judged (fixed before running), with `diagnose_chains.py`:** do the chains agree; are the remaining sampled parameters believable (Tier 2 and log prior); how large sigma_m becomes (a larger value means more misfit admitted openly); the canopy swing (lowest / highest LAI); and the average misfit by season (script 05, part 5). Preferred: the least-fixed version whose remaining parameters are believable. All versions are reported.

### D55. Expansion test: a cold limit on photosynthesis, "DALEC2cold" (Decided as a test, 2026-10-16)
**Why.** F18 and F20: with model error the computation is fine, but every version of DALEC2 makes up for the missing winter shutdown of photosynthesis through whatever canopy parameters are free, and fixing them all only worsens the fit. The predictive checks show the cause directly: too much uptake in spring, too much release in winter. Gelman et al. (2020, Sections 2.2, 7 and 10) treat adding a missing process as the normal response to a misfit that matters for the questions asked.

**What.** Stettz et al. (2022), DALEC2 at an evergreen subalpine conifer forest: GPP is multiplied by

    g(Tmin) = 0 at or below t_shutdown;  1 at or above t_full;  (Tmin - t_shutdown) / (t_full - t_shutdown) in between

with both thresholds calibrated. Here Tmin is the daily minimum (the model is daily; Stettz et al. used monthly values).
- **Priors** (`COLD_PRIORS` in `priors.py`, revised 2026-10-16 to follow the official code): `t_shutdown` lies in -15 to 0 degC and `t_full` in 0 to 15 degC, the ranges of `PARS_INFO_831.c`. The ranges do not overlap, so `t_full` can never fall below `t_shutdown`. The official prior is flat over each range; here each threshold is a Beta(2, 2) fraction of its range (same range, smooth, most weight in the middle, no hard walls for NUTS; D33). (First version, replaced: `t_shutdown ~ Normal(-7, 3)`, `t_full = t_shutdown + gap`, gap about 7 to 24 degC.)
- **Code:** `cold_multiplier` (`acm.py`), `run(..., cold=...)` (numpy), `build_forward(..., cold=...)` (PyTensor), `build_model(..., cold_limit=True)`, `--cold` in script 05; the posterior checks and `diagnose_chains.py` apply it too. **Tested:** PyTensor equals numpy to 1e-10; the gradient with respect to a threshold matches a numerical slope; the full log posterior equals the hand computation; the ramp is right at both thresholds (58 tests in total).

**First test:** a quick fit (real NEE 1997-1998, model error, nothing fixed), judged like F18 and F20: do the chains agree; is the forest believable (Tier 2, log prior, peak LAI, `ceff`, `f_fol`); the model error sigma_m; the seasonal misfit. **If the canopy becomes believable**, the cold-limit model goes through the full workflow before it is used for results: prior predictive check (script 03), fake-data check (script 04), then the full real-data fit. If not, it is reported as tried.

**Source checked: the paper and the official code (2026-10-16).** Stettz et al. (2022), Eqs. (1) and (2), and their code release CARDAMOM-FluxVal v1.0 (Yang et al. 2021, Zenodo doi:10.5281/zenodo.4904195), which contains CARDAMOM v2.2. The cold limit is in `C/projects/CARDAMOM_MODELS/DALEC/DALEC_831/DALEC_831.c`, lines 166 to 188, marked "SGS 7/17/19: Temp scaling factor" (S. G. Stettz). Checked line by line:
- **Same as ours:** the ramp (0 below `Tminmin`, 1 above `Tminmax`, straight line between); the temperature is the time step's minimum (`MET[m+1]`; monthly in their study, daily here); the factor multiplies GPP only (line 188).
- **Priors:** `PARS_INFO_831.c`, parameters 33 and 34: 258.15 to 273.15 K and 273.15 to 288.15 K, i.e. -15 to 0 and 0 to 15 degC, flat. Adopted as ranges (above).
- **A quirk not copied:** the code stores the thresholds as whole numbers (`int Tminmin = pars[33] - 273.15;`), which rounds them to whole degrees towards zero. Harmless for their Metropolis sampler, which uses no gradients; NUTS needs smooth gradients, so our thresholds stay continuous.
- **Not copied either:** DALEC_831 is built on a richer DALEC than DALEC_400 (a soil-water bucket limiting GPP, evapotranspiration, fire, a sampled labile lifespan, respiration relative to the site mean temperature with a precipitation factor). Only the cold limit is taken; our base model stays DALEC_400 (D21 to D31).

**Not part of DALEC2 as published;** a stated model change, reported next to the unchanged model (F17 to F20).

### D56. Is the thin canopy forced by the data, or only allowed by them? (Decided as a test, 2026-10-16)
**Why.** F21: the cold limit fixed the winter shrinking of the canopy and improved the fit, but the canopy stays thin (peak LAI about 1) with very efficient leaves (`ceff` about 110) and heavy leaves (`lma` about 435). Two readings are possible. (a) NEE barely tells a thin, efficient canopy from a full, ordinary one (equifinality), and the fit drifts to the thin one for small gains; then holding the canopy at believable values costs little fit. (b) The data clearly prefer the thin canopy, because something else in the model is missing or wrong; then holding the canopy costs a lot of fit.

**What.** Two quick fits (real NEE 1997-1998, model error, cold limit), compared with the free cold fit of F21 (log-likelihood about -509):
- `--cold --fix lma,c_lf,ceff`: canopy held at believable values (S2 of D54 plus the cold limit).
- `--cold --fix ceff`: only leaf efficiency held at the REFLEX median 10.

**How it is judged (fixed before running):** the drop in log-likelihood against the free cold fit, together with sigma_m.
- Drop under about 10: reading (a). NEE cannot pin down the canopy; the believable canopy is kept (fixed, or with priors that carry the site knowledge), and this is an RQ3 result.
- Drop over about 30: reading (b). The data prefer the thin canopy; the next suspects are the canopy part of the model (the ACM response to LAI and ceff), or adding LAI observations as a second data stream.
- In between: no clear answer from quick fits; decided on full fits with a proper comparison (PSIS-LOO, Vehtari et al. 2017).
(For scale: S2 without the cold limit lost about 110 against S1, F20.) This is a diagnostic comparison of quick fits, not a formal model choice.

### D57. Priors of the cold-limit model: `ceff` returns to its first prior, D38 (Tested 2026-10-16/17: FAILED, so not adopted; see F24)
**Why.** F22: with `ceff` near 10, ACM photosynthesis hardly rises with sunlight, and the fit escapes to an implausible canopy. D42 had lowered `ceff` from D38's median 16.5 (95% 7.4 to 37, built from p11 x needle N before any data were seen) to median 10 (95% 5 to 20). Its only stated reason: with no cold limit, higher values make annual GPP implausibly high, because the canopy photosynthesises through the winter. The cold-limit model removes that reason, so the narrowing has to be checked again for that model.

**What.** In the cold-limit model only, `ceff ~ LogNormal(median 16.5, sdlog 0.41)` (D38). Everything else as before. The main model keeps D42.
- Code: `COLD_MODEL_PRIOR_CHANGES` and `priors_for(cold)` in `priors.py`; `build_model(..., cold_limit=True)` uses them. Each fit now saves its prior list with the results (`priors` note in the .nc file), so the parameter table and `diagnose_chains.py` always use the priors the fit was made with. Fits saved earlier carry no note; they all used the main priors.

**Check before any fit (prior predictive, script 03 `--cold`, no tower NEE used).** Same 1,000-draw design and D40 checks as F13 to F15, with the soft foliage constraint (D44) as fitted. **Rule, fixed before running:** use the D38 `ceff` prior in the cold model if Tier 2 is at least 50% and Tier 1 at most 7% (the main model was accepted at 6.4%, F15; with about 800 effective draws the Monte Carlo error of such a rate is about 0.9 percentage points). Otherwise keep D42 for the cold model too (option B of F22). The script also runs the same draws with the D42 `ceff` prior, for comparison.

**Outcome (F24): FAIL on Tier 1 (18.7%, mostly daily GPP above 20). Under the rule, the cold-limit model keeps the D42 `ceff` prior**; `COLD_MODEL_PRIOR_CHANGES` is now empty. The code that saves each fit's prior list stays.

**Planned next step if it had passed:** a quick fit of the cold model (`--cold`), judged like F21 and F22. If it is believable, the cold model goes through the fake-data check and the full fit (D55).

**Gelman et al. (2020) on this step:**
- Section 2.4: prior predictive checks "provide a way to refine the model without using the data multiple times". That is why the decision rests on the prior predictive check and on the D38 and D42 sources, not on the fit.
- Section 7.3: "Priors must be specified for each model within a workflow." The cold-limit model is a new model, so its priors are set and checked on their own.
- Section 7.3 also notes that "often as the model becomes more complex, all of the priors need to become tighter." Here one prior is loosened. The reason is specific: the tight `ceff` prior did the job of the missing cold limit, and the cold limit now does that job itself.
- Section 6.3: sensitivity analysis by "refitting the model with multiple priors". The cold fit with the D42 prior (F21) is kept and reported as that comparison.
- **Stated openly:** this change was considered after seeing that the fit prefers high `ceff` (F21, F22). It is not a free choice of the best-fitting prior: it returns to a prior written down before any data were seen, for a reason recorded in D42 at the time, and it must pass a check that uses no tower data.

### D58. Diagnostic: estimate the ACM light-use efficiency `a7` (Decided as a test, 2026-10-17)
**Why.** F24 suggests that the thin canopy compensates for an ACM light-use efficiency that is too high for this site: `a7` = 7.19 is fixed at the value fitted at the temperate Loobos site (D22), and shrinking the canopy is the model's only way to lower it. If that is right, letting `a7` be estimated should give a normal canopy with a lower `a7` and an equally good fit.

**What.** Cold-limit model (D55) with the main priors (D57 outcome), plus `a7` sampled: `acm_a7 ~ LogNormal(median 7.19298, sdlog 0.5)`, 95% about 2.7 to 19 (`ACM_PRIORS` in `priors.py`). No boreal value exists to centre it on, so it is centred on the fixed value and kept wide; the F24 calculation (about 45% of 7.19, i.e. 3.2) lies inside it. Options `--free-a7` in scripts 03 and 05 (with `--cold`); results named `..._cold_a7...`.
- Code: `acm.gpp(..., a7=)`, `run(..., a7=)`, `build_forward(..., a7=)`, `build_model(..., free_a7=True)`; tests check numpy against PyTensor, the gradient with respect to `a7`, and the full log posterior against a hand computation.
- `diagnose_chains.py` shows `a7`; its "log prior" leaves out the `a7` and cold-threshold priors, as before for the thresholds, so it stays comparable with F17 to F24.

**Steps and rules, fixed before running:**
1. **Prior predictive check** (`03 --cold --free-a7`), same rule as D57: Tier 2 at least 50% and Tier 1 at most 7%, with the soft foliage constraint. If it fails, no fit.
2. **Quick fit** (`05 --cold --free-a7 --quick`), compared with the free cold fit of F21 (log-likelihood -508.7, peak LAI 1.07, `ceff` 110):
   - **Supports the hypothesis** if all hold: the chains agree; peak LAI at least 2; `ceff` inside the official CARDAMOM range 5 to 50; `a7` median below 5; log-likelihood not more than 10 below F21.
   - **Contradicts it** if the canopy stays thin (peak LAI below 1.5) with `a7` near or above 7.
   - **Not identifiable** if the chains disagree along the `a7`-canopy trade-off, or the `a7` interval is about as wide as its prior: then NEE cannot tell light-use efficiency from canopy size, which is reported as an RQ3 result.
   - Anything else: reported as unclear.
3. **A diagnostic, not a model choice.** Even if supported, estimating `a7` is not automatically adopted. Adoption would be a separate decision and would need the full workflow (fake-data check with `a7` sampled, full fit, comparison with the cold model by PSIS-LOO).

### D59. Fake-data check of the cold-limit model (Decided, 2026-10-17)
**Why.** The cold-limit model (D55) with the main priors passed its prior predictive check (F24). Before its full fit can be trusted computationally it needs its own fake-data check (Gelman et al. 2020, Section 4.1), like the main model (D48, F16).

**What.** Script 04 with `--cold`: the "true" forest is the first fully plausible prior draw of the cold model, including two true thresholds; fake NEE is made with the cold limit; the fitted model has the cold limit. Fake noise is RANDUNC only, so the fit uses the RANDUNC-only observation model, as in D48. Sampler: dense mass matrix with target_accept 0.95 (D53). Results named `fake_data_check_cold`.

**Criteria (fixed before running):** the four D49 criteria. The truth criterion keeps the same allowance of misses as D49 (4): with the two thresholds there are 26 quantities, so at least 22 of 26 must contain the truth (with 90% intervals, a perfect fit misses 5 or more of 26 in about 11% of runs, against about 9% for 5 or more of 24, so the test is almost equally strict).

### D60. Figures and tables from a saved fit, for talks and the thesis (Decided, 2026-10-19)
**What.** `scripts/talk_figures.py <label>` reads `results/<label>.nc` and makes, without sampling again: a summary table (csv, plus two booktabs halves for slides), trace plots (R), rank plots (S), prior against posterior (T) and posterior correlations (U) for every sampled quantity; with `--predict` (real-data fits) it runs 200 posterior draws forward for the NEE runs (V), the pools (W) and the sources of uncertainty (X). `scripts/show_data_columns.py` prints the data columns the model uses, for the data slide; it saves nothing.

**Choices.**
- **Rank plots next to trace plots** (Vehtari et al. 2021): a trace plot can look fine while chains mix slowly; rank plots show it as uneven bars.
- **Prior against posterior uses the priors stored in the fit** (D57), so a fit made with other priors is drawn against its own.
- **Correlations on the sampler's scale** (log or logit, D32). A strong correlation means NEE fixes a combination of the two quantities, not each one (RQ3).
- **Figure X (RQ2) is descriptive.** Parameter uncertainty is the spread of the model's daily NEE across posterior draws, averaged by month over all days; measurement error is the mean RANDUNC on the days used; model error is the median `sigma_m`. The three are standard deviations: they do not add up (the total predictive spread is the square root of the sum of their squares).
- **Figure W:** pools every 365 days from 1 January 1997, so "start of each year" drifts by one day per leap year.
- **Slide table:** mean, 5%, 95%, R-hat and ESS, the smaller of bulk and tail as in the D49 rule; values failing the D49 rule are bold. Rates below 0.01 are shown in units of a power of ten (for example 10^-5 d^-1). The csv keeps every column.
- **Layout:** sized for a 16:9 slide (drawn at about twice the slide size, with large fonts) and without titles: the slide or the thesis caption says what a figure shows. Chain colours were checked for colour-blind safety as a set.

### D61. Testing on the held-back years and comparing models: rules fixed before reading the full fits (Decided, 2026-10-21; closes OPEN-1)
**When written.** After the two full fits of 2026-10-19 finished (jobs 2369912, cold model; 2369913, main model with the dense setting) and before any of their results beyond the run time were read. Fixing the rules first stops the choice of test days or of the margin from being tuned to the answer.

**Which fits.** Only fits that pass the D49 computation checks (no divergences, R-hat below 1.01, ESS at least 400 bulk and tail) are tested or compared. The `a7` model (D58, F28) is not part of the comparison: it has not been through its own prior and fake-data checks.

**Main test: the held-back years 2011 to 2014 (Gelman et al. 2020, Section 6.2; the plan of D4).** For each fit, 1,000 posterior draws are run forward from 1 January 1997 straight through 31 December 2014 with the drivers of every day, so the pools at the start of 2011 come from the fitted years and the model has never seen the test years. Measured on the test days:
1. **Log predictive density (the decision score):** for each test day, log of the average over draws of Normal(observed NEE | model NEE, sqrt(RANDUNC^2 + sigma_m^2)), summed over the test days. It rewards a model that is both close to the data and honest about its uncertainty.
2. Coverage of the 90% predictive interval, overall and by season.
3. RMSE and mean bias (observed minus model) of the posterior mean.
4. NEE summed per test year over the test days (like for like).

**Test days (closes OPEN-1).** Primary: the same rule as for fitting (D12): NEE present, QC at least 0.75 and RANDUNC present, 810 days in 2011-2014. Every measure uses these days. Secondary: RMSE and bias also on the 1,418 days that pass QC alone (RANDUNC is not needed for those two), reported alongside, not used for the decision.

**When is one model clearly better.** The difference in the decision score, with the standard error from the day-by-day differences (as in Vehtari et al. 2017). Neighbouring days are correlated (F11), so the standard error is widened by sqrt((1 + r) / (1 - r)), with r the correlation between the differences on consecutive test days (only when r is positive). A model is clearly better when its score is higher by more than twice this widened standard error. Otherwise there is no clear difference: both are reported, and the preferred model is the one whose forest is believable (share of draws passing the D40 Tier 2 checks over 1997-2010, and peak LAI against the site's 2.5 to 4).

**Also reported, not used for the decision:** PSIS-LOO on the 4,413 assimilable days of 1997-2010 (Vehtari et al. 2017), from the same 1,000 draws, with the Pareto k diagnostic. It is stated as optimistic: leaving out one day while its correlated neighbours stay in makes prediction look easier than it is.

**Code:** `scripts/06_evaluate.py <fit> [<fit> ...]`; settings in the `evaluation` section of `config/default.yaml`. Nothing is sampled again.

### D62. The main model is not fitted again; the cold-limit model becomes the base model (Decided, 2026-10-23)
**Why.** The main model (DALEC2 without the cold limit) has failed D49 on real NEE twice: with the diagonal setting (F27) and with the dense setting (F30). The second failure is not a tuning problem: the chains sit in two separate modes, and **no chain is a believable forest** (Tier 2 0%, peak LAI about 1, log prior far below every prior draw). Gelman et al. (2020, Section 5.1, the "folk theorem") read this as a problem of the model, not of the sampler, and Section 5.9 warns against changing sampler settings until the problem goes away. Smaller steps (a higher `target_accept`) help with divergences, not with separate modes.

**What.**
1. **The main model is not run again.** F27 and F30 are kept and reported as results: on 14 years of real NEE, DALEC2 with the ACM light-use efficiency fixed and no cold limit cannot be fitted to a believable forest at FI-Hyy.
2. **The D61 comparison "main against cold" is reported as not possible.** D61 compares only fits that pass D49, and the main model's fit does not. This is stated in the thesis, not hidden.
3. **The cold-limit model (D55) is the base model from now on, called M2** (26 quantities plus `sigma_m`). It passed its prior predictive check (F24) and its fake-data check (F26). Its full fit failed only on 7 divergences (F31), so it is fitted again with `target_accept` 0.99, the planned response of D49 and F31 (job "ta99": `fit_model_error_cold_dense_ta99`).
4. **The question "which model" moves to the reduced cold models of D63.**

### D63. Reduced cold models: holding fixed what the fake data showed NEE cannot learn (Decided, 2026-10-23, before any of these runs)
**Why.** The supervisor asked to keep constant the parameters that, in the fake-data check, got almost no information from NEE. Gelman et al. (2020) support this: Section 5.2 (simplify the model when the computation is hard), Section 4.1 (the fake-data check shows which quantities the data can and cannot learn), and Section 7 (building a model up step by step and comparing the steps). A parameter that NEE cannot learn adds width to the sampler's search without adding anything the data can check, and it can create the long curved valleys where divergences happen.

**Which quantities.** Chosen only from the cold model's fake-data check (F26), using the width ratio (posterior width divided by prior width; 1 means NEE taught nothing). The real-data fits were not used to choose. The ladder:

| Model | Held fixed | Free quantities | Reason for this step |
|---|---|---|---|
| M2 | nothing (the cold model, D62) | 26 | the base |
| MA | starting wood pool `c_woo_0` | 25 | the least informed (0.87); also the only "not informed" quantity in the main model's check (F16) |
| MB | MA + starting soil pool `c_som_0` | 24 | the other large, slow pool (0.75): 14 years of NEE barely see where these two pools started |
| MC | MB + soil and wood turnover `theta_som`, `theta_woo` | 22 | the rates of the same two slow pools (0.75 each) |
| M0 | MC + the four growth shares (together), `h`, `cr_onset`, `lma`, `f_auto` | 14 | everything not "strongly informed" (width ratio 0.5 or more, the D49 boundary) |

Notes on the ladder:
- **Ordered by pool, not strictly by width ratio.** The share sent to the roots (`f_roo`, 0.86) is about as uninformed as `c_woo_0`, but the four growth shares come from one Dirichlet prior and sum to 1, so one share cannot be fixed alone. They are fixed together, at the Dirichlet mean, and only in M0.
- `c_roo_0` (0.495) and `theta_roo` (0.487) are just under 0.5, so they stay free in M0. The 0.5 boundary was set in D49, before F26, and is not moved.
- **Each model is nested in M2:** it is M2 with some quantities held at one value. So the ladder asks how much the fit and the predictions lose (or gain) when those quantities are no longer learned.

**Fixed values.** Each fixed quantity is held at its prior median. The four growth shares are held at the mean of their Dirichlet prior, times (1 - `f_auto`) as in D28; when `f_auto` is fixed it is also at its prior median. The minimum turnover (`theta_min`) follows from the fixed values as before. The priors of the fixed quantities drop out of the posterior; nothing else in the model changes.

**Fake-data check of each reduced model (Section 4.1, as D48 and D59).** The true forest is drawn from the reduced prior, so the fixed quantities take their fixed values in the truth too: the check asks whether the sampler and the remaining quantities work, not whether fixing is right (that is asked on real data, below). Same design and the same D49 criteria. The number of misses allowed is set by one rule for every model size: the smallest m such that, if each truth falls outside its 90% interval with probability 0.1, more than m misses happen in at most 12% of runs (`fake_data.miss_tail`). This gives 4 misses for 22 to 26 quantities, the same as D49 (20 of 24) and D59 (22 of 26), and 3 for the 14 of M0 (at least 11 of 14).

**Sampler for all ladder runs:** dense mass matrix, `target_accept` 0.99 (the setting M2 now needs, D62), given on the command line, so file names end in `_dense_ta99`. The fake-data checks use the RANDUNC-only error model, as in D48; the real fits use the main error model (L1).

**How the model is chosen (fixed now, before any result).**
1. A model takes part only if its fake-data check passes and its real-data fit passes D49. A model that fails is reported with its failure.
2. The models that take part are scored on the held-back years with the D61 decision score (log predictive density on the 810 test days).
3. Every model whose score is within 2 widened standard errors of the best score (D61) counts as equally good.
4. **Among the equally good models, the one with the fewest free quantities is chosen** (simplest first).
5. **Safeguard:** a model whose forest is not believable (share of draws passing the D40 Tier 2 checks over 1997-2010 below 50%) is not chosen over an equally good believable one, even if it is smaller.
6. Every model on the ladder is reported (scores, D49 results, believability), not only the chosen one. PSIS-LOO on the fitted years is reported as in D61, not used for the choice.

**Not done, stated as limitations.**
- No separate prior predictive check for the reduced models. Fixing quantities at their medians only removes spread from a prior that already passed (F24), but the reduced prior is not exactly a part of the full one, so this is an assumption. It can be checked cheaply later if needed.
- The real-data posteriors may show other quantities as uninformed. They are not used to build a second ladder: that would choose the model from the same data it is tested against.

**Code.** `fix_sets` in `config/default.yaml`; option `--fix-set MA` (MB, MC, M0) in scripts 04 and 05; fixing in `src/dalec/priors.py` (`fixed_values`, `apply_fixed`) and `src/dalec/inference.py` (the shares become a fixed value instead of a Dirichlet); checks in `tests/test_fix_sets.py`, including that the log density of M0 equals the hand computation of the full model minus the priors of the fixed quantities.

### D64. Literature priors and the simplest model that represents the forest (Decided 2026-10-27, before any literature value was read; the supervisor's suggestion, OPEN-9)
**Why.** The supervisor asked for the simplest model that can represent the forest, with posterior values from the literature used as priors, especially for quantities NEE does not inform (such as the starting wood pool). This follows Gelman et al. (2020): simplify (Section 5.2), build up step by step (Section 7), and use outside information in the priors and check how much the priors matter (Sections 6.3 and 7.3). Famiglietti et al. (2021) found that a model of intermediate complexity predicted best at FI-Hyy and five other flux sites.

**Source (OPEN-9 rules 1 and 2).** Bloom et al. (2016), CARDAMOM 2001-2010 (doi:10.7488/ds/316): the same DALEC2 code (`DALEC_400`), fitted to MODIS LAI and soil carbon (HWSD), with no NEE, so its data are separate from ours. The 1 degree grid cell that contains FI-Hyy (61.85 N, 24.29 E): 61 to 62 N, 24 to 25 E. Famiglietti et al. (2021) is not used as a prior, because it fitted FI-Hyy NEE from about 1999 to 2003, which overlaps our data.

**Which quantities.** Only the weakly informed ones (the M0 set: width ratio 0.5 or more in the cold fake-data check, F26) that CARDAMOM reports: `f_auto`, the four growth shares, `theta_woo`, `theta_som`, `lma`, `c_woo_0` and `c_som_0`. CARDAMOM does not report `h` and `cr_onset`; they keep the site priors. All strongly informed quantities keep the site priors (D34 to D39): NEE decides them anyway.

**How a CARDAMOM number becomes a prior (fixed now):**
1. **Centre:** the cell's median (or exp(log-mean) where only log statistics are given).
2. **Spread:** the cell's spread on the same scale (log sd for positive quantities, sd for fractions), **doubled**. This counts the CARDAMOM information as a quarter of its weight (a "power prior" with weight 1/4), because it describes a 1 degree cell from satellite LAI and a global soil map, not this stand.
3. **Rates from residence times:** `theta_woo` = 1 / (365.25 x wood residence time in years). `theta_som` = 1 / (365.25 x soil residence time x F), where F is the mean over 1997-2010 of exp(Theta T_mid) at the prior median Theta = 0.0336. Our `theta_som` is the rate at 0 degC (D20, D35), while CARDAMOM's residence time includes the effect of temperature. If the units turn out not to be years, the units are converted; the rule stays.
4. **Pools:** `c_woo_0` and `c_som_0` from the 2001-2010 mean pools (an approximation: these slow pools change little between 1997 and 2010).
5. **Shapes:** log-normal for rates, pools and `lma`; beta for `f_auto` (mean = the cell median, sd doubled); Dirichlet for the four shares (means = each share of GPP divided by the sum of the four; concentration set so that the average sd of the four shares equals the doubled CARDAMOM sd, each divided by 1 - `f_auto`).
6. The values go into `src/dalec/priors.py` (`PRIOR_SETS["lit"]`) next to the CARDAMOM numbers they came from, and are printed and recorded here before any run. They are not changed after any check or fit.
7. **Details settled after seeing the file layout, before printing any value at the cell (2026-10-27):** the cell is the one centred at 61.5 N, 24.5 E (its box is 61 to 62 N, 24 to 25 E). The files give mean, median, sd and the 5th, 25th, 75th and 95th percentiles; the residence-time files give a "log-based" mean and sd. For LMA and the pools, which have no log statistics, the log sd comes from the 25th and 75th percentiles (the quartiles: the most robust pair). Whether the residence-time "log-based sd" is the sd of the log or its exponential is judged over all land cells against the quartiles of the same file, not at our cell. Units: pools in t C ha-1, times 100 for g C m-2; residence times in years; LMA in g C m-2. Code: `scripts/07_literature_priors.py`.

**The ladder of literature models ("L"), mirroring D63.** Cold model, dense mass matrix, target_accept 0.99. A fixed quantity is held at the median of its literature prior; the shares at the literature Dirichlet mean; `h` and `cr_onset` (in L0) at their site medians, as in M0.

| Model | Held fixed | Free | Partner in D63 |
|---|---|---|---|
| L2 | nothing (literature priors on the 10 quantities above) | 26 | M2 |
| LA | `c_woo_0` | 25 | MA |
| LB | + `c_som_0` | 24 | MB |
| LC | + `theta_som`, `theta_woo` | 22 | MC |
| L0 | the M0 set | 14 | M0 |

**Steps for each (as D63):** a prior predictive check of the literature priors (script 03, the D40 numbers; **reported, not a gate**, because changing the literature values after the check would be tuning them); a fake-data check with the D63 misses rule; a full fit (D49); the test years (D61).

**How the model is chosen ("the least complexity that represents the forest"), across every eligible model of D63 and D64:**
1. Eligible: the fake-data check passes and the full fit passes D49.
2. **Represents the forest:** at least 50% of posterior draws pass the D40 Tier 2 checks over 1997-2010.
3. Among those, the one with the fewest free quantities; a tie goes to the higher test score.
4. Its test score is compared with MA by the D61 rule. If it is not clearly worse, it becomes the main model. If it is clearly worse, both are reported: MA as the best predictor of NEE, the chosen one as the simplest believable forest, used for the pools and fluxes of RQ2.
5. If no model represents the forest, MA stays the main model, and D64 is reported as "literature priors did not give a believable forest", with how far they moved it.
Every L model is also compared with its D63 partner (site values against literature values).

**Expected before running (written so the result cannot be re-told afterwards):** the weakly informed quantities hardly change NEE, so their new priors should change the fit little. The exceptions may be `lma` and the foliage share, which shape the canopy: values from an LAI-constrained analysis may pull the canopy towards a believable one, at some cost in fit, because NEE prefers the thin canopy (F22).

**Cost:** 10 jobs on Roihu (5 fake-data checks, 5 fits), each up to about a day on 4 cores, plus one short prior check: roughly 3,000 to 6,000 billing units of the 57,000 left.

### D65. What to do with D64 after reading the literature values (Proposed 2026-10-27; to agree with the supervisor; no fit has been run)
**Why a new decision.** F35 shows that the literature values, built by the D64 rules, are centred far from the stand, carry little information after widening, and include a Dirichlet with all weights below 1, which breaks D33. Running D64 as written is still possible, but not running it is a change of plan, so it is recorded here with its reasons, before any fit.

**Options:**
- **A. Run D64 as written** (10 jobs, about 1 day). Expected: most prior draws impossible, the fake-data checks may not even find a plausible true forest (0.7% of draws in the dry run), and NUTS may struggle with the Dirichlet. Likely uninterpretable.
- **B. Run only the prior predictive check of the literature priors on the site weather** (minutes), report F35 and that check, and stop D64 there (recommended). The thesis then states that the only independent DALEC2 posterior for this location describes a 1 degree, deciduous-like, peat-soil ecosystem, so site measurements are the better source of priors for this stand; and it uses F35's agreement (thin canopy, high Rh) in the Discussion.
- **C. In addition to B, an exploratory test of the supervisor's own example:** hold the starting wood pool (MA) and then the wood and soil pools (MB) at the CARDAMOM medians (10,760 and 23,610 g C m-2), with every other prior at its site value (prior sets `lit_wood`, `lit_wood_soil`). 4 jobs (2 fake-data checks, 2 fits), about a day. It asks whether the value a weakly informed quantity is fixed at matters. Labelled exploratory, because it was chosen after reading the values. Expected: the wood pool changes little; a soil pool 4 times the site value forces the soil turnover or Rh to change.
**Recommended: B now, C if the supervisor wants a run.**

### D66. The litter-to-soil rate gets its prior directly from the published range; the share `h` is removed (Decided by the author, 2026-10-31; replaces the `h` part of D35)
**What.** `theta_min` (litter to soil organic matter, d-1, at 0 degC) is sampled directly with a log-normal prior whose 95% range is the published range of Bloom and Williams (2015, Table 1), 1e-5 to 1e-2 d-1: median 3.16e-4 (the geometric middle), sdlog 1.762. This is the D35 fallback rule for a parameter without a site value; the shape is smooth, like every other prior (D33). The share `h` (D35, from RothC) is no longer used anywhere. In the M0 fix set, `theta_min` takes the place of `h`.

**Why.** The author prefers a prior taken from the DALEC2 literature over one built from another model (RothC), which would need its own defence in the thesis. The model equations are unchanged; only this prior changes. The old prior via `h` had a similar centre (median 2.8e-4) but a range about 10 times narrower.

**Consequence: every result that used the old prior is run again,** with the same settings and rules (D49, D61, D63): the prior predictive check of the cold model (D40/D57 numbers), the five fake-data checks (M2, MA, MB, MC, M0) and the five full fits, then the test years (D61) and the D63 choice. F26, F31, F32 and F34 stay in this file as the record of the runs with the `h` prior; the thesis reports the new runs. The fix sets are not re-chosen: MA to MC do not involve this quantity, and M0 keeps its role ("everything not strongly informed").

### D67. The allocation prior written as three beta steps (Decided by the author, 2026-10-31; replaces the form of D28 and D34, not their values)
**What.** The four shares of GPP that go to the living pools are built step by step, in the order of the official DALEC2 code (D28): foliage takes a fraction `frac_fol` of the growth carbon (GPP - Ra), labile a fraction `frac_lab` of what is left, fine roots a fraction `frac_roo` of what is left after that, and wood the rest. Each fraction has a beta prior: `frac_fol ~ Beta(13.06, 72.59)`, `frac_lab ~ Beta(13.06, 59.53)`, `frac_roo ~ Beta(15.26, 44.27)`. The sampler samples these three fractions; the four shares are computed from them.

**Why.** Easier to explain, and it follows the official code's own split. **It is exactly the same prior as before** (D34): each step is a beta of its own weight against the sum of the weights still to come (checked with 400,000 draws: every share's 2.5%, 50% and 97.5% points and the correlations between shares agree to the third decimal). The sources and the reasoning of D34 are unchanged.

**Consequences.**
- The sampled quantities are now 23 for DALEC2 and **25 for the cold model, the same count as the model's parameters** (17 process parameters, 6 starting pools, 2 thresholds): `f_auto` and the three fractions stand for the four allocation parameters of Bloom and Williams (2015, Table 1), and `theta_min` is sampled directly (D66). The fake-data checks judge these 25; the misses rule (D63) allows 4, so at least 21 of 25 for the full cold model.
- The M0 fix set holds the three fractions (instead of "the four shares together"), so M0 still has 14 free quantities.
- The literature prior set (D64, on hold under D65) is written in the same step form; its values are the exact step equivalent of the earlier ones.
- Together with D66, every run of the D63 ladder is repeated with the new code (jobs submitted on 2026-10-31 with the D66 code were cancelled before they finished).

### D68. A second true forest, judged together with the first, and a longer, more careful sampler (Decided 2026-11-01, after F36 and before any new run)
**Why.** In F36 the misses of the fake-data checks come from the true forest itself: it sits in the tails of the priors (LMA 121, below the prior's 5% point 123; roots fraction 0.357, above the prior's 95% point 0.353; `ceff` 6.15 near its lower end), and the misses are near misses towards the prior's centre (`f_auto` missed by 0.0004). A 90% interval covers the truth 90% of the time on average over truths drawn from the prior, so one truth in the tails says little. D48 already allowed "a second truth ... if the first shows problems". Separately, `theta_min` mixes slowly (moderate correlations with `theta_lit` and `c_lit_0`, both -0.46), and MB and MC had one divergence each.

**What (fixed now):**
1. **Second truth:** the second fully plausible prior draw of the same search (`--truth 2`), not a chosen one.
2. **Truth criterion judged over both truths pooled:** all free quantities of both runs count together, and the misses allowed follow the D63 rule for the pooled number (M2: at least 42 of 50; MA: 41 of 48; MB: 39 of 46; MC: 35 of 42). Each run must still pass the computation checks on its own. Code: `scripts/08_pool_fake_checks.py`.
3. **Sampler for every run of this round (both fake-data checks and the fits):** 2,000 kept draws per chain (was 1,000), `target_accept` 0.995 (was 0.99), dense mass matrix, warm-up unchanged (1,000). File names end in `_dense_ta995_d2000`. This changes the computation only, not the model or the rules.
4. **Models:** M2, MA, MB and MC. M0 is not run again: it failed its fake-data check clearly in both rounds (F32, F36) and is reported as such.
5. **Then** D61 (test years) and D63 (choice) as written, on the models that pass.
**If a model still fails,** it is reported as failing; the rules are not changed again in this round.

### D69. Prior sensitivity of the chosen model by power-scaling (Decided by the author, 2026-11-03, before any result)
**Why.** RQ1 asks how much NEE constrains each quantity. The width ratio (posterior width / prior width) shows how much the data narrowed a prior, but not whether the prior still steers the answer. Gelman et al. (2020, Section 6.3) recommend checking how sensitive the posterior is to the prior. Refitting with other priors costs about 27 h per fit; power-scaling needs no new sampling.

**What.** For the model chosen under D63 only: the prior and the likelihood are each raised to a power near 1 (ArviZ `psense_summary`, default alphas 0.99 and 1.01), and the saved posterior draws are re-weighted (Kallioinen et al. 2024). Per quantity this gives a prior sensitivity and a likelihood sensitivity, read with the ArviZ threshold 0.05:
- prior high, likelihood low: the prior decides (the data say little);
- both high: prior-data conflict;
- both low: the data decide.

**What is scaled.** The prior of every sampled quantity, `sigma_m` included. The soft foliage constraint (D44) is not scaled: it is not a prior on one quantity, so it stays the same in every perturbation. The likelihood is the NEE likelihood (L1) on the 4,413 assimilable days.

**Reported, not a gate.** The result describes the chosen model; it does not change the choice or the priors (changing priors after seeing it would tune them to the data). A strong prior-data conflict is reported and discussed (for example, `theta_min` near the low end of its range, F36).

**Addition (2026-11-04, after reading the joint result, F40; stated openly):** the joint table scales all priors at once, so a quantity can show prior sensitivity only because it is correlated with another quantity whose prior moves. A second table scales each quantity's own prior on its own (`--by-quantity`, ArviZ `prior_var_names`), from the saved log prior and log likelihood; no model run. It is a reading aid for the same result, reported next to the joint table, and like it changes no prior and no choice.

**Code.** `scripts/11_prior_sensitivity.py <fit>` (rebuilds the model from the fit name, adds the log prior and log likelihood of each saved draw, saves `results/<fit>_psense.nc` and `figures/<fit>_prior_sensitivity.csv`); checks in `tests/test_prior_sensitivity.py`, including that the added terms equal the model's own prior and likelihood at a draw.

## Part 3. Likelihood

### L1. Observation model: what sigma is (Decided, 2026-10-08; main model changed 2026-10-16)

**Change on 2026-10-16: the main model is now option C (RANDUNC plus model error).**

    NEE_t ~ Normal( model_t , sqrt(RANDUNC_t^2 + sigma_m^2) ),   sigma_m ~ HalfNormal(1)

**Why it changed.** On the real NEE, option A left the sampler with separate modes and pushed the parameters to physically impossible forests (F17). Option C, planned below as the sensitivity run before any data were seen, removed the separate modes and brought the parameters much closer to the priors (F18). The fitted model error, about 0.6 g C m-2 d-1, is larger than RANDUNC on most days: DALEC2's daily misfit is larger than the measurement error, so treating RANDUNC as the whole error claimed a precision the model cannot have. This follows Gelman et al. (2020, Sections 5.1, 7.1 and 10.5: an added error term "to allow for model misspecification"). **Option A becomes the sensitivity run** and is reported (F17), not hidden. The fake-data check (D48) still uses option A, because its fake noise is RANDUNC only.

**The original main model (option A, 2026-10-08), kept for the record: sigma is RANDUNC.**

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

- Option C: sigma_t = sqrt( RANDUNC_t^2 + sigma_m^2 ), with sigma_m ~ HalfNormal(1) estimated jointly. Adds one shared model-error term. (Main model since 2026-10-16, see above.)
- Support for option C if needed: Gelman et al. (2020, Section 10.5), golf putting example (known per-point variance plus an added independent error term "to allow for model misspecification"); Fer et al. (2018), observation-error parameters estimated during calibration.
- Compare the two on: posterior widths, winter share of the likelihood weight, R-hat, divergences and number of posterior modes.

Not chosen:

- A constant sigma (option B): ignores that sigma grows with flux size (R^2 about 0.47 in the previous analysis).
- sigma = k * RANDUNC: rescales every day equally, so the relative weights between days stay the same (the reason Gelman et al. give for rejecting the beta-binomial).

Implementation (planned): the observation model will be a config setting (`randunc_only` for the main fit; `randunc_plus_model_error` for the sensitivity run), so switching needs no code change. It is not in `config/default.yaml` yet; it is added together with the PyMC model.

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

**F5. Missing RANDUNC.** 1,138 of 6,940 days (16.4%) in the full file have no RANDUNC, including all of 2014. Within the calibration block, 365 days lack it, spread over several gaps (shown in Figure A, which also prints every gap of 7 days or more), not one whole year.

**F6. Evaluation block.** 608 of 1,461 evaluation days (41.6%) are excluded only because RANDUNC is missing: all of 2014 and the first half of 2011.

**F7. Daily temperature extremes are not in the daily file (resolved by D18).** DALEC2's photosynthesis (ACM) uses the daily temperature range, but the FLUXNET daily file has no daily maximum or minimum. `TA_F_DAY` and `TA_F_NIGHT` are only stand-ins; the previous repository found they understate the range by about 4.8 C and derived true extremes from the half-hourly file. This is the one model input the daily file cannot supply directly.

**F8. Driver completeness.** In the full record, the model drivers are complete except CO2 (0.3% empty: the 19 days at the start of 1996, the reason 1996 is excluded, D4). Inside the calibration block every model input is complete (F10). Soil temperature and soil water are 1 to 5% empty; net radiation 21%, diffuse shortwave 46%, soil heat flux 46%, diffuse photosynthetic light 77%.

**F9. Temperature-related cautions from earlier DALEC users (literature, 2026-10-12).**
- Bloom and Williams (2015) chose test sites "with a ≤ 3 months of recorded below-freezing soil temperatures" and little water stress, because "these criteria reflect the current capabilities of DALEC2". FI-Hyy, a boreal site, has a longer frozen season (the previous repository's notes say about four to six months of soil frost; to be confirmed from the soil temperature data). Using DALEC2 here is a deliberate scope choice and a stated limitation.
- Stettz et al. (2022), DALEC2 at a subalpine evergreen forest (US-NR1), report that baseline DALEC2 has no cold-temperature control on spring photosynthesis: GPP is "limited primarily by incoming shortwave radiation". They added a scalar on daily minimum temperature (GPP gradually reduced below about 6.0 C, fully stopped below about -7.1 C), which improved spring GPP (r = 0.88 vs 0.47). This matches the early spring onset found in the previous repository and supports showing spring in Figure I.

**F10. First real-data run of the driver figures (2026-10-12).**
- The half-hourly safety check (D18) passed: every day has 48 valid half-hours, and the half-hourly mean matches the daily `TA_F` within 0.05 C on every day.
- True daily temperature range averages 6.3 C over 1997-2010. All model inputs are complete in the calibration block (no missing days).
- Spring onset with the 10-day rule (D16): median day 100, range 66 to 137 over the calibration years; evaluation years 93 to 111. See OPEN-3.

**F11. Days are not independent (Figure K).** NEE anomalies are correlated with the next day (+0.37), two days later (+0.22) and a week later (+0.08); they reach zero by about two weeks. If days were independent, values would stay within +/-0.03. The Gaussian likelihood (L2) assumes independent days, so it treats the data as more informative than it is and the posterior may be too narrow. Recorded as a limitation; an error model with day-to-day correlation (for example AR(1), where each day's error carries part of the previous day's) is a possible later sensitivity analysis.

**F12. Prior predictive check, round 1 (real data, 2026-10-13; D40).** 1,000 draws, seed 20261013.
- **Tier 1 FAIL:** 34.7% of draws physically impossible (target at most 5%; the previous repository's published-range priors: 84%). Reasons: daily GPP above 20 (34.1%), annual GPP above 3,000 (20.8%), LAI above 12 (14.1%).
- **Tier 2 FAIL:** 29.7% pass every check (target at least 50%). Per check: GPP 47.7%, Rh 63.1%, NEE 63.1%, peak LAI 54.8%, foliage ratio 48.5%, wood ratio 64.3%, soil ratio 55.2%.
- **Priors linked to failure:** high `ceff` (+1.03), low `lma` (-0.86), low `c_lf` (-0.66), high `f_lab` and `f_fol`, low `f_auto`. Together: the **canopy feedback** (more leaf area, more GPP, more allocation to leaves). Respiration, soil and timing priors show little link.
- **Plausible draws, medians:** GPP 1,560 (site about 1,030), Rh 488 (site 290 to 370), NEE -248, peak LAI 3.7 (site about 2.5 to 4), foliage end/start 0.83, wood 1.49, soil 1.04.
- **Frost switch on (D24):** Tier 1 29.5%, Tier 2 37.8%; better, still failing.
- **NEE seasonal cycle (Figure L, observed for context only):** winter and autumn match the observed cycle closely; summer observed values sit inside the middle half; **in spring nearly all prior draws turn to net uptake by early March, about six weeks before the observed cycle (mid-April)**. Because this appears across the whole prior, it is structural (no cold limit on photosynthesis; F9, Stettz et al. 2022), not something a prior can fix.
- **Next (D40, "if the check fails"):** revisit only the canopy priors (`ceff` first, then `lma`, `c_lf`, allocation to leaves), each with a stated source. **Done:** `ceff` in D42, `c_lf` in D43, then the soft constraint D44; `lma` and the allocation priors were left unchanged (D45).

**F13. Prior predictive check, round 2 (only `ceff` changed, D42; real data, 2026-10-13).**
- **Tier 1:** 9.3% impossible (round 1: 34.7%; target at most 5%), still FAIL. Daily GPP above 20: 7.7%; LAI above 12: 6.2%; annual GPP above 3,000: 3.2%.
- **Tier 2:** 42.6% pass every check (round 1: 29.7%; target at least 50%), still FAIL. The weakest check is now the **foliage pool, end / start (55.7%)**, with many draws near 0: the canopy **collapses** in a fifth to a quarter of draws.
- **Plausible draws, medians:** GPP 1,061 (site about 1,030: now matches), Rh 425 (site 290 to 370), NEE -68, peak LAI 3.1, foliage ratio 0.58, wood 1.28, soil 1.00.
- **Links:** Tier 1 (runaway canopy): low `lma`, high `ceff`, low `c_lf`, low `f_auto`, high `f_lab` and `f_fol`. Tier 2 (mainly canopy collapse): high `c_lf` (+0.40) first.
- **Reading:** the canopy is unstable in both directions. Needles are lost at rate `c_lf` and replaced by allocation `f_fol + f_lab` times GPP, which itself depends on leaf area (`C_fol / lma`) and `ceff`. These four priors are set independently, so some combinations replace needles faster than they fall (runaway) and others slower (collapse): the joint effect Gelman et al. (2020, Section 7.3) warn about. `c_lf` appears in both failure modes.
- **Frost switch on:** Tier 1 6.9%, Tier 2 42.2%.
- **Re-run with the full parameter table (same seed): identical numbers**, so the check is reproducible. Full list adds `f_woo` (+0.43, Tier 1; low `f_auto` leaves more carbon for every plant pool) and, for Tier 2, `theta_woo` (+0.22) and `c_som_0` (-0.21). **Caution on reading links:** with about 93 failing draws, a standardised difference has a standard error of roughly 0.11, and 24 quantities are compared at once; links below about 0.3 (here `c_som_0`, `d_fall`, `theta_woo`) may be chance and are not acted on.
- **Figure L:** winter and autumn still match; the spring early uptake remains (structural, F12); the prior median summer uptake (about -1.3) is now weaker than observed (about -3), but the observed values stay inside the 5 to 95% band.

**F14. Prior predictive check, round 3 (only `c_lf` changed, D43; real data, 2026-10-13).**
| | Round 1 | Round 2 | **Round 3** | Target |
|---|---|---|---|---|
| Tier 1, impossible | 34.7% | 9.3% | **8.6%** | at most 5% |
| Tier 2, pass all | 29.7% | 42.6% | **47.6%** | at least 50% |
| foliage end/start inside band | 48.5% | 55.7% | **61.8%** | |
| GPP median (plausible draws) | 1,560 | 1,061 | **1,122** | site about 1,030 |
- With the frost switch on: Tier 1 5.8%, Tier 2 46.0%.
- Medians: Rh 433, NEE -93, peak LAI 3.4, foliage 0.63, wood 1.31, soil 1.00.
- Links: Tier 1 still the canopy group (`ceff` +1.15, `lma` -1.11, `f_auto` -0.81, `f_lab` +0.78, `f_fol` +0.68, `c_lf` -0.67); Tier 2 links all below 0.3 (`c_lf` +0.29), i.e. no single prior stands out any more.
- **Reading:** both targets are now narrowly missed, and narrowing single priors has reached diminishing returns (Tier 1: 34.7, 9.3, 8.6). The remaining failures are canopy runaway and collapse, which is exactly what Bloom and Williams (2015) remove with their ecological dynamic constraints EDC 7 (no pool grows faster than 10% per year) and EDC 8 (no pool decays with a half-life under 3 years). This project dropped EDCs because hard accept/reject rules break NUTS (carried-over decision). **What followed:** the official code turned out to use a different, steady-state form, and its smooth version was added for the foliage pool (D44, F15).

**F15. Prior predictive check with the soft foliage constraint (D44; real data, 2026-10-13).** Same 1,000 draws as round 3, re-weighted by the constraint (effective number of draws 808; 1.0% of draws get almost no weight).
| | Round 3 (independent priors) | **With D44** | Target |
|---|---|---|---|
| Tier 1, impossible | 8.6% | **6.4%** | at most 5% (narrow miss) |
| Tier 2, pass all | 47.6% | **63.1%** | at least 50% (PASS) |
| foliage end/start in band | 61.8% | 81.3% | |
| GPP in band | 72.5% | 83.7% | |
| peak LAI in band | 74.3% | 85.0% | |
- Remaining Tier 1 failures: daily GPP above 20 (5.8%), LAI above 12 (3.1%), annual GPP above 3,000 (2.1%).
- With the frost switch on as well: Tier 1 4.8%, Tier 2 64.7% (both pass); reported alongside only (D24, D40), not used to choose.
- Part of the foliage improvement is expected by construction (D44 caveat); GPP and LAI checks improved too, which the constraint does not target directly.

**F16. Fake-data check, full run (D48, D49; Roihu job 2278889, 2026-10-14/15).** Fake NEE from a known "true" forest on the real 1997-2010 weather, gaps and RANDUNC; 4 chains x (1,000 warm-up + 1,000 kept).
- **D49 criteria: all PASS.** No divergences (0 of 4,000); largest R-hat 1.002; smallest effective sample size 1,963 (bulk) and 2,304 (tail); truth inside the 90% interval for **20 of 24**, exactly the minimum set in advance.
- **The four misses:** `f_fol` (only just outside: truth 0.0484, interval from 0.0487), `theta_roo` (truth 0.0016, interval 0.0018 to 0.0040), `c_roo_0` (truth 372, interval 177 to 346) and `c_som_0` (truth 4,439, interval 5,057 to 9,160). Four or more misses happen by chance in about 21% of perfect fits, so the count itself is not alarming. **But the misses are not random:** they come in pool and rate pairs. For roots, the fit chose a smaller pool turning over faster than the truth; for soil, a larger pool (with a slightly slower rate, still inside its interval). NEE sees mainly the flux leaving a pool (about rate x pool), so "small pool, fast rate" and "large pool, slow rate" look alike. The implied fluxes still differ from the truth by about 18% (roots) and 29% (soil), so other parameters share the compensation. Reading: these pools and rates are **only weakly identifiable from NEE alone**, even when the model is perfect. The soil truth also sat in the low tail of its prior (about -0.8 on Figure N's scale), and with weak information the posterior stays nearer the prior centre. This is a result for RQ1, not a code error.
- **How much NEE narrowed each prior** (posterior width / prior width): strongly informed (below 0.5): `f_auto`, `theta_lit`, temperature exponent (0.03), `d_onset` (0.08), `cr_onset`, `d_fall`, `cr_fall`, `c_lf`, `ceff` (0.07), `c_lab_0`, `c_fol_0`, `c_lit_0`. Partly informed: the four allocation shares, `theta_roo`, `theta_woo`, `theta_som`, `h`, `lma`, `c_roo_0`, `c_som_0`. Not informed: `c_woo_0` (0.92). Caveat: one truth only (D48); a second truth would show whether this pattern is general.
- **Cost:** 29.4 h on 4 CPUs (52.9 s per draw). The step size settled at about 0.005 and almost every draw used the maximum 1,023 gradient evaluations (92 of 4,000 hit the maximum tree depth, which wastes time but does not bias the result). Likely cause: the posterior is very narrow in some directions (temperature exponent, `ceff`, `d_onset`) and has long sloping ridges in others (the pool and rate pairs), which a diagonal mass matrix cannot straighten out. The first job with a 24 h limit (2278575) timed out, so D50's limit was raised to 48 h.
- **What it allows:** the code and sampler recover a known answer, so a real-data fit with the same settings can be trusted computationally (Gelman et al. 2020, Section 4.1).

**F17. First contact with the real NEE: the quick fit leaves the believable region (D51; Roihu job 2316333 and `diagnose_chains.py`, 2026-10-16).** Two short chains on 1997-1998 (622 assimilable days); a quick run, so nothing is final, but the pattern is clear.
- **The chains disagree completely** (largest R-hat 2.54, effective sample size about 3), unlike the fake-data quick run (R-hat 1.03).
- **Neither chain is a believable forest.** A typical draw from the prior has a log prior of about -22, and only 1 in 1,000 is below -37. Chain 0 sits at -384 and chain 1 at -998.
  - Chain 0 passes Tier 1 but not Tier 2: needles living about 1.3 years (`c_lf` 0.79; prior 3 to 5 years), `ceff` 247 (prior 5 to 20), a starting foliage pool of 10 g C m-2 (prior about 590), peak LAI 1.2 (site 2.5 to 4), temperature exponent 0.097 (Q10 about 2.6; prior 1.2 to 2.0), starting soil pool 1,313 (prior about 6,100).
  - Chain 1 is physically impossible (peak LAI 32, annual GPP 2,868).
- **Which chain the model prefers:** chain 0, by 1,899 log units; it also fits the data better, by 1,114.
- **Reading.** Two problems, one on top of the other:
  1. **A model problem (Gelman et al. 2020, Section 5.1).** Both chains started from believable prior draws (D45) and moved away from them, and the better chain paid about 360 log units of prior to gain fit. So the data prefer an unbelievable forest over any believable one. The likelihood is worth thousands of log units and the priors hundreds, because RANDUNC is small and the model cannot match the daily NEE within it (structural misfits F12, F13). The data overwhelm the priors.
  2. **Separate modes (Section 5.6).** Chain 1 sits in a region about 1,900 log units worse and did not leave it: a computational symptom of the same rugged posterior.
- **The fake-data check (F16) is what makes this reading possible:** the same code and sampler recovered a known forest when the model was true, so the failure here points to the mismatch between DALEC2 and the real data, not to a bug.
- **Correction:** the first version of `diagnose_chains.py` called chain 0 "believable" because it passed Tier 1. That was too weak a test; the script now also requires Tier 2 and a log prior within the range of real prior draws.
- **Next:** open, see OPEN-8.

**F18. Quick fits of the OPEN-8 options (real NEE 1997-1998, 2 short chains each; `diagnose_chains.py`, 2026-10-16).** Log prior of a typical prior draw about -22 (1 in 1,000 below -37).

| Quick run | Chains agree? | Log prior (both chains) | Model error sigma_m |
|---|---|---|---|
| RANDUNC only (F17) | no, 1,899 apart | -384 / -998 | (none) |
| + frost switch | no, 903 apart | -340 / -1,040 | (none) |
| + model error | yes, 2.6 apart | -89 / -90 | 0.59 g C m-2 d-1 |
| + model error + frost | yes, 0.8 apart | -76 / -78 | 0.57 |

- **Model error removes the separate modes** and brings the parameters much closer to the priors: the error model was the cause of the computational problem (Gelman et al. 2020, Section 5.1). The fitted model error, about 0.6 g C m-2 d-1, is larger than RANDUNC on most days, i.e. DALEC2's daily misfit is bigger than the measurement error.
- **The frost switch alone changes nothing;** with model error it helps a little.
- **What remains is the same in both model-error runs:** annual GPP about 1,020 (inside the site range 952 to 1,104) produced by an implausible canopy: peak LAI about 1 (site 2.5 to 4), `ceff` about 110 (prior 5 to 20), `lma` 400 to 460, needles living 2.1 to 2.6 years (prior 3 to 5), `f_auto` 0.34 (prior 0.37 to 0.65), Q10 about 2.0 to 2.3.
- **The canopy shrinks each winter.** Lowest / highest leaf area within a year: 0.40 (model error) and 0.50 (with frost), against 0.69 to 0.83 (median 0.77) for believable prior forests (computed on synthetic weather; this ratio depends mainly on the leaf timing parameters). **Interpretation (a hypothesis, consistent with the numbers):** DALEC2 has no cold limit on photosynthesis (F12; Stettz et al. 2022), so the fit imitates the winter shutdown of the real forest by making the canopy partly deciduous: small, short-lived needles, very efficient leaves. Total GPP right, canopy wrong: the "right answer for the wrong reasons" of RQ3.

**F19. Speed test with a dense mass matrix (D52; Roihu job 2317355, fake data, 2026-10-15).**
- **Time: 4.8 h against 29.4 h** (about 80 instead of about 1,000 gradient evaluations per draw; step size 0.049 against 0.005). D52 criterion 2: PASS.
- **Same answer:** every posterior median within 2% of the diagonal run's 90% interval width (largest `f_woo`, 1.9%); interval widths within 7 to 9%; the same four truths outside their intervals. D52 criterion 3: PASS.
- **But 13 divergences in 4,000 draws** (0.3%), so D49 criterion 1 fails and D52 criterion 1 fails. R-hat (largest 1.005) and effective sample sizes (at least 592) pass.
- **Verdict under the rule fixed in advance: dense is not adopted; the diagonal setting stays.** The divergences most likely come from the much larger step size in a few curved regions; the usual remedy is a higher `target_accept` (smaller steps). That is a new test, D53, not a reinterpretation of D52.

**F20. Quick fits of the simpler models S1, S2, S4 (D54; real NEE 1997-1998, model error on; Roihu jobs 2350529 to 2350531, 2026-10-16).**

| Version | Chains agree? | Log prior | Peak LAI | Annual GPP | `ceff` | `f_fol` | sigma_m | Log-likelihood |
|---|---|---|---|---|---|---|---|---|
| nothing fixed (F18) | yes | -89 | 0.9 | 1,020 | about 110 | 0.08 | 0.59 | -589 |
| S1: `lma`, `c_lf` fixed | yes | -87 | 1.0 | 1,180 | about 105 | **0.015** | 0.65 | -643 |
| S2: + `ceff` = 10 | yes | -43 | **8** | **1,840** | 10 (fixed) | 0.07 | 0.80 | -755 |
| S4: S1 + three pools fixed | yes | -84 | 1.0 | 1,170 | about 104 | 0.016 | 0.65 | -643 |

(Log prior of a typical prior draw about -22; fixed parameters enter it at their median. Site values: peak LAI 2.5 to 4, GPP 952 to 1,104; `f_fol` prior median 0.07.)
- **No version is a believable forest; the compensation moves to whatever is still free.** With `lma` and `c_lf` fixed, the fit keeps the canopy thin by starving the foliage of allocation (`f_fol` 0.015, a fifth of the prior) and keeps `ceff` about ten times its prior. With `ceff` also fixed, the canopy instead grows to LAI 8 and GPP to about 1,800, and the fit to the data gets clearly worse (log-likelihood 110 lower, model error 0.80).
- **Fixing the poorly informed pools (S4) changes nothing:** the problem is in the canopy, not the pools.
- **The same seasonal misfit in all three** (predictive checks, part 5): coverage of the 90% interval about 89% overall, but 100% in winter and about 79% in summer, so one constant model error is too wide in winter and too narrow in summer. Average observed minus model: **spring +0.17 to +0.29** (the model takes up too much CO2 in spring: the D51 expectation, confirmed) and **winter -0.19 to -0.29** (the model releases too much in winter).
- **Caveat on the canopy-swing number** (F18): in these short fits the foliage pool also grows strongly over the two years (soft foliage term -1.4 to -3.4), which lowers "lowest / highest LAI in a year" on its own. `diagnose_chains.py` now also prints the foliage pool's end / start ratio, so trend and seasonality can be told apart.
- **Reading.** Simplifying removes freedom but not the cause: whichever canopy parameters stay free absorb the missing cold shutdown, and fixing all of them turns the gap into a worse fit instead. This points to a structural gap in DALEC2 at this site (the model cannot shut photosynthesis down in winter and early spring), which is what the expansion route (a cold limit on GPP, Stettz et al. 2022) would address. It is also an RQ3 result in its own right: NEE alone cannot tell canopy size, leaf efficiency and allocation apart when the model structure is wrong.

---

**F21. Quick fit with the cold limit (D55; real NEE 1997-1998, model error on, nothing fixed; 2026-10-16).** Compared with the same fit without the cold limit (F18):

| | Without cold limit (F18) | With cold limit |
|---|---|---|
| Chains agree? | yes | yes (log posterior 1.1 apart) |
| Log prior | -89 / -90 | -68 / -68 |
| Model error sigma_m | 0.59 | 0.51 |
| Lowest / highest LAI in a year | 0.40 | 0.61 (believable prior forests 0.69 to 0.83) |
| `f_fol` | 0.08 (F20; starved to 0.015 in S1 and S4) | 0.071 (prior median about 0.07) |
| Peak LAI | about 1 | 1.07 |
| `ceff` | about 110 | 110 |
| `lma` | 400 to 460 | 425 to 445 |
| Annual GPP | about 1,020 | about 1,090 |

Thresholds: no GPP below about -6.6 degC, full GPP above about 1.2 degC (daily minimum).
- **What the cold limit fixed:** the fit is better (sigma_m 14% smaller), the canopy no longer shrinks much in winter, foliage allocation is at its prior (it was starved in S1 and S4), and the log prior improves by about 21.
- **What it did not fix:** the canopy is still thin, with very efficient and heavy leaves. Still not a believable forest (Tier 2: 0%; log prior far below the 1-in-1,000 prior floor of -37).
- **Correction to F18:** F18 read the whole implausible canopy as imitating the missing winter shutdown. That explains the winter shrinking, but not the thin, efficient canopy, which stays with the cold limit in place. The cause of the thin canopy is open; D56 tests it.
- Also seen: Q10 about 2.3 (`temperature_exponent` 0.084, prior median 0.034, i.e. Q10 1.4); `c_lf` about 0.35 (needles about 2.9 years, prior 3 to 5).
- **Seasonal misfit (average observed minus model, g C m-2 d-1; F20 values without the cold limit in brackets):** winter DJF -0.12 (-0.19 to -0.29), spring MAM -0.04 (+0.17 to +0.29), summer JJA +0.06, autumn SON -0.01. **The spring misfit is gone; the winter misfit is smaller but still there.** This matches what D51 expected from a missing cold limit.
- **Coverage of the 90% predictive interval:** all 89.9%, DJF 98.5%, MAM 90.3%, JJA 80.7%, SON 91.6%. The interval is too wide in winter and too narrow in summer. One shared model error (sigma_m) cannot be right in both seasons: the model's misfit is larger in summer, when fluxes are large. A model error that grows with the size of the flux is a candidate for later (random flux errors grow with flux size; Richardson et al. 2006), tried only after D56, one change at a time.

**F22. Result of D56: the data strongly prefer the thin canopy (real NEE 1997-1998, model error and cold limit on; Roihu jobs 2354407 and 2354408, 2026-10-16).**

| Quick run | Log-likelihood | Drop against free cold fit | sigma_m | Peak LAI | Annual GPP | `ceff` | `lma` | Tier 1 impossible |
|---|---|---|---|---|---|---|---|---|
| Free (F21) | -508.7 | | 0.51 | 1.07 | 1,090 | 110 | 435 | 0% |
| `ceff` held at 10 | -650.7 | 142 | 0.66 | 14.0 | 2,095 | 10 | 87 | 96% |
| `lma`, `c_lf`, `ceff` held | -659.1 | 150 | 0.67 | 11.2 | 1,918 | 10 | 185 | 48% |

(Means of the two chains; the chains agree within 2 units in both runs.)
- **Verdict under the D56 rule: reading (b).** Both drops are far over 30. NEE does not merely allow the thin, efficient canopy; it clearly prefers it.
- **With ordinary leaves (`ceff` 10) the fit does not settle on a normal forest either.** It builds a huge canopy (LAI 11 to 14) and about twice the site's GPP (site 952 to 1,104), then removes the extra carbon with high respiration (`f_auto` 0.60 to 0.65). So with this model there are two ways to fit NEE: a thin canopy of very efficient leaves with the right GPP, or a huge canopy of ordinary leaves with double GPP. Neither is the real forest.
- **The official CARDAMOM ranges confirm that the free fit is outside what DALEC users accept** (`PARS_INFO_831.c`, same values as for the base model): `ceff` 5 to 50 (free fit 110), leaf mass per area 5 to 200 g C m-2 (free fit 435). Temperature factor 0.018 to 0.08 (Q10 up to 2.2); all fits here sit at its upper edge (0.075 to 0.084).
- **Why, from the ACM itself (computed with `acm.gpp`, a summer day at Hyytiala: Tmax 20, Tmin 10 degC, CO2 370, day 180).** GPP on a dull day (8 MJ m-2 d-1) and a bright day (28 MJ m-2 d-1):

| Canopy | Dull day | Bright day | Bright / dull |
|---|---|---|---|
| LAI 3, `ceff` 10 (our prior centre) | 5.2 | 6.7 | 1.30 |
| LAI 3, `ceff` 30 (needle-N value, D42) | 8.7 | 14.1 | 1.63 |
| LAI 3, `ceff` 50 (official upper limit) | 9.9 | 17.6 | 1.78 |
| LAI 11, `ceff` 10 (held fit) | 10.2 | 16.3 | 1.60 |
| LAI 1.07, `ceff` 110 (free fit) | 5.3 | 11.7 | 2.19 |

  With `ceff` near 10, ACM photosynthesis is capped by its leaf-efficiency (carboxylation) term, so GPP hardly rises with sunlight. NEE evidently varies with sunlight much more than that, and the fit has only two ways to get that response: very efficient leaves on a thin canopy, or a huge canopy. **The thin canopy is the model's way of making GPP follow light.** Hyytiala's own GPP model is a light-use-efficiency model (PRELES; Peltoniemi et al. 2015), which fits this picture. Still to confirm on the data: residuals of the held fit against daily radiation (Gelman et al. 2020, Section 6.1).
- **Link to D42.** D42 lowered the `ceff` prior from D38's median 16.5 (7.4 to 37) to median 10 (5 to 20), and recorded the tension: needle nitrogen points to about 30. Its stated reason was that higher values make annual GPP implausibly high **in a model with no cold limit** (the canopy photosynthesises all winter). The cold-limit model removes that reason, so for the cold model the D42 narrowing has to be checked again by its own prior predictive check (proposed as D57). This is said openly: the posterior is now known to prefer high `ceff`, so any change must rest on the D38 and D42 sources and on the prior predictive check, not on the fit.

**F23. Second speed test: dense mass matrix, target_accept 0.95 (D53; Roihu job 2344869, fake data, 2026-10-16).**
- **D49 criteria: all PASS.** 0 divergences in 4,000 draws, 0 at the maximum tree depth; largest R-hat 1.004 (`h`); smallest effective sample size 1,038 (bulk) and 1,219 (tail); truth inside the 90% interval for 21 of 24 (outside: `theta_roo`, `c_roo_0`, `c_som_0`; `f_fol`, just outside in F16, is now inside).
- **Time: 6.2 h against 29.4 h** for the diagonal run (about 119 instead of about 1,000 gradient evaluations per draw; step size about 0.036). Under the 20.5 h limit.
- **Same answer:** compared with the diagonal run (F16) on the prior scale, every posterior median moved by at most 2.2% of the diagonal run's 90% interval width (largest `f_roo`), far under the 10% limit; interval widths 0.96 to 1.07 times the diagonal ones.
- **Verdict under the rule fixed in advance: adopted** for full runs (D53). Raising target_accept from 0.9 (F19, 13 divergences) to 0.95 removed the divergences, as expected, at a cost of about 1.4 h.

**F24. Prior predictive check of the cold-limit model with the D38 `ceff` prior (D57; Roihu interactive job 2358821, real weather, no NEE used, 2026-10-17).** 1,000 draws; the same draws re-run with the D42 `ceff` prior.

| With the soft foliage constraint (as fitted) | D38 `ceff` (median 16.5) | D42 `ceff` (median 10) | Target |
|---|---|---|---|
| Tier 1, impossible | **18.7%** (daily GPP above 20: 18.5%) | 4.0% | at most 7% (D57) |
| Tier 2, plausible | 60.7% | 63.6% | at least 50% |
| Median annual GPP of plausible draws (unweighted) | 1,096 | 652 | site 952 to 1,104 |
| Median peak LAI of plausible draws (unweighted) | 3.2 | 2.45 | site 2.5 to 4 |

- **Verdict under the rule fixed in advance: FAIL.** The cold-limit model keeps the D42 prior (D57 option B). The strongest link to failure is high `ceff` (+1.24), then low `lma` (-0.88).
- **What fails is the daily peak, not the yearly total.** With D38 the plausible draws have about the right annual GPP and LAI, but almost one draw in five has summer days above 20 g C m-2 d-1 (observed daily GPP at Hyytiala stays far below that).
- **The cold model with the D42 prior passes** (Tier 1 4.0%, Tier 2 63.6%), slightly better than the main model (F15: 6.4%, 63.1%). So the cold-limit model with unchanged priors has a sound prior predictive check; its first quick fit is F21.

**What F22 and F24 together point to: ACM's fixed light-use efficiency.** In the ACM, light-limited photosynthesis is `e0 x radiation`, with `e0 = a7 LAI^2 / (LAI^2 + a9)` and `a7` = 7.19, `a9` = 2.10 fixed (fitted at the temperate Loobos site, D22). `e0` is 2.5 at LAI 1.07 and 5.8 at LAI 3. A summer day computed with `acm.gpp` (Tmax 20, Tmin 10 degC; dull 8, bright 28 MJ m-2 d-1):

| Canopy | Dull | Bright | Bright / dull |
|---|---|---|---|
| Free cold fit: LAI 1.07, `ceff` 110 | 5.3 | 11.7 | 2.19 |
| LAI 3, `ceff` 110, `a7` as fixed | 10.9 | 21.0 | 1.93 |
| LAI 3, `ceff` 16.5, `a7` as fixed | 6.8 | 9.8 | 1.43 |
| LAI 3, `ceff` 30, `a7` x 0.45 | 5.2 | 10.8 | 2.06 |

- The data ask for GPP that rises strongly with light (F22) but stays moderate on bright days. With `a7` fixed, a normal canopy can give one or the other: efficient leaves on a normal canopy overshoot (21 on a bright day, as in F24), ordinary leaves barely respond to light. **Shrinking the canopy is the only way the model has to lower its light-use efficiency.**
- The last row shows the alternative: a normal canopy (LAI 3) with `ceff` 30, the value needle nitrogen suggests (D42), and `a7` at about 45% of the Loobos value gives the same response as the free fit.
- **Hypothesis (not yet tested on data):** the thin canopy compensates for an ACM light-use efficiency that is too high for this boreal site. This also explains why the fit could not be repaired by priors on the canopy (F20, F22) or on `ceff` (F24). A possible test: let `a7` be estimated (D58, proposed).

**F25. Prior predictive check of the cold-limit model with `a7` sampled (D58 step 1; Roihu interactive job 2358821, real weather, no NEE used, 2026-10-17).** Same 1,000-draw design; the same draws re-run with `a7` fixed.
- **With the soft foliage constraint (as fitted): Tier 1 4.3%, Tier 2 62.0%. D58 rule: PASS**, so the quick fit goes ahead. With `a7` fixed, the same draws give 4.0% and 63.6%: sampling `a7` hardly changes what the priors imply about the forest.
- Without the constraint Tier 2 is 37.2% (39.4% with `a7` fixed); the constraint mainly removes draws whose foliage pool runs down, as in F15.
- Tier 1 failures link most to high `ceff` (+1.33), low `lma` (-1.13) and high `a7` (+0.85): the same "too much light-use" direction as F24.
- Also submitted on 2026-10-17: the fake-data check of the cold-limit model (D59), Roihu job 2359179.

**F26. Fake-data check of the cold-limit model: PASS (D59; Roihu job 2359179, dense mass matrix, target_accept 0.95, 17 h 6 min, finished 2026-10-18).**
- **D49 criteria: all pass.** 0 divergences of 4,000 draws, no draw at the maximum tree depth; largest R-hat 1.0051 (`f_lab`); smallest ESS bulk 2,280, tail 2,420; truth inside the 90% interval for **24 of 26** (rule: at least 22).
- **The two misses.** `cr_onset`: truth 106 d, interval 40 to 82 (the same parameter is only partly informed in F16 too: NEE sees the timing of leaf growth more than its speed). `cold_t_full`: truth 5.57 degC, interval 4.84 to 5.41, so missed by 0.16 degC with a very narrow interval (posterior width 4.5% of the prior's). Two misses are about what chance gives with 26 quantities at 90%; the narrow `cold_t_full` interval is noted as something to watch on real data (many days near the threshold make it look better known than it may be).
- **Both cold thresholds are strongly informed** (`t_shutdown` 13.9 to 13.2 below zero, width 14% of the prior's). NEE can identify the cold limit, which supports estimating it (D55).
- **Informed by NEE:** 14 strongly, 12 partly, 0 not informed (main model F16: the starting wood pool was not informed).
- **Cost:** about 500 gradient evaluations per draw (step size 0.0094), 17 h against 6.2 h for the main model with the same settings (F23). A full real-data fit of the cold model should be planned for well over a day.
- **Consequence (D55, D59):** the cold-limit model is cleared, computationally, for its full fit on real NEE.

**F27. Full fit of the main model on real NEE: FAILS the D49 computation criteria (D51; Roihu job 2344870, 29 h 11 min, finished 2026-10-17). Results not interpreted (D51 step 1).**
- Largest R-hat 1.081 (`f_lab`); smallest ESS bulk 40, tail 84 (rules: below 1.01, at least 400). About 1,001 gradient evaluations per draw with a step size of 0.0011: the sampler was close to its maximum tree depth on most draws, the same sign of hard geometry as the first fake-data check (F16, D52). The log confirms the cause: it ran with the **diagonal** mass matrix and target_accept 0.9 (submitted before D53 made the dense matrix the default), had **29 divergences**, and **3,845 of 4,000 draws hit the maximum tree depth**. All three D49 computation checks fail.
- **Not interpreted, recorded only:** the summaries look like the quick fits (peak LAI about 1.0, annual GPP about 1,050, Rh 508 above the site range 290 to 370, coverage of the 90% interval 89% overall with 99.9% in winter and 78% in summer). With R-hat 1.08 these numbers may change in a converged run.

**F28. D58 quick fit: `a7` estimated (real NEE 1997-1998, model error and cold limit; Roihu job 2359287, 2 chains x 100 draws, 2026-10-17). Verdict under the D58 rules: UNCLEAR (4 of 5 conditions met; `ceff` just above its limit). Partly supports the hypothesis.**
- **Each D58 condition:** chains agree, roughly (largest R-hat 1.03, ESS 75, no divergences; usual for a quick fit). Peak LAI median **2.23** (90%: 1.65 to 2.92; rule: at least 2): **met** by the median. `ceff` median **54** (40 to 76; rule: inside 5 to 50): **not met**, just above. `a7` median **3.16** (2.57 to 4.08; rule: below 5): **met**. Log-likelihood **-514.7 and -514.4** per chain (`diagnose_chains.py`), 5.8 below F21 (-508.7; rule: not more than 10 below): **met**. Per chain `ceff` is 64.6 and 54.6, both above 50. One condition fails, so by the rule fixed in advance the result is "unclear", not "supports". This is reported as it stands; the rule is not changed after seeing the result.
- **What moved:** compared with the free cold fit (F21), peak LAI went from 1.07 to 2.23 and `ceff` from 110 to 54. `a7` came out at 44% of its fixed value 7.19, close to the F24 calculation (about 45%, i.e. 3.2) made before the fit. `a7` is strongly informed (posterior width 28% of the prior's), so NEE can identify it: the "not identifiable" outcome did not happen.
- **What else moved:** the cold thresholds changed a lot (`t_shutdown` about -6 instead of -12.5 degC, `t_full` about 1.3 instead of 6 degC) and the temperature exponent rose to 0.083 (Q10 about 2.3). Heterotrophic respiration (460) is still above the site range (290 to 370). So freeing `a7` changes how the model splits the work between light use, cold and decay, not only the canopy.
- **Still far from the priors:** `diagnose_chains.py` calls both chains "not believable": log prior about -53, while only 1 prior draw in 1,000 is below -36.6 (this log prior leaves out the `a7` and cold-threshold priors). Besides `ceff` (55 to 65, prior 5 to 20), needles live about 2.7 years (`c_lf` 0.37; prior 3 to 5 years), `f_auto` is 0.31 (prior 0.37 to 0.65) and the temperature exponent is about 2.4 times its prior median. Tier 2 plausibility is high (92% and 96%) because those checks look at the fluxes and LAI, not at the parameters.
- **Reading:** the fixed `a7` explains the thin canopy (canopy and `a7` both move as predicted, and NEE identifies `a7`), but not the rest of the tension between the data and the priors. The fit to the data does not improve (5.8 worse); freeing `a7` moves the solution to a believable canopy at a small cost in fit.
- **Seasonal fit:** 90.5% of days inside the 90% interval; winter -0.12, spring -0.01, summer +0.03, autumn -0.01 (similar to the cold fit, F21).
- **Status:** a diagnostic only (D58 step 3). Even a clear "supports" would not adopt `a7` without the full workflow.
- **Submitted on 2026-10-19 (72 h limit):** full fit of the cold-limit model on real NEE (Roihu job 2369912, `fit_model_error_cold`), and the main model again with the dense mass matrix and target_accept 0.95 (job 2369913, `fit_model_error_dense_ta95`), so that both models can be compared once converged.

**F29. Full fit of the cold-limit model on real NEE: narrowly FAILS D49; results not interpreted (Roihu job 2369912, 32.1 h, finished 2026-10-20).**
- 0 divergences, but largest R-hat 1.0123 (`theta_som`) and smallest ESS 337 bulk (`lma`); rules: below 1.01, at least 400. 3,969 of 4,000 draws hit the maximum tree depth (about 1,023 gradient evaluations per draw, step size 0.0018).
- **Cause found:** the job ran with the **diagonal** mass matrix and target_accept 0.9, not the D53 setting. The D53 change to `config/default.yaml` had not reached the repository on Roihu when the job was submitted (the `git pull` of 2026-10-21 brought it). The main-model rerun (F30) was not affected because its settings were given on the command line.
- Not interpreted (D51 step 1). Recorded only: the summaries resemble the cold quick fit (peak LAI about 1.05, annual Rh 555 above the site range).

**F30. Main model on real NEE with the D53 setting: FAILS D49 badly; results not interpreted (Roihu job 2369913, dense mass matrix, target_accept 0.95, 33.7 h, finished 2026-10-20).**
- 17 divergences, largest R-hat **2.21** (`d_fall`; also `c_lf` 1.97, `cr_fall` 2.05, `cr_onset` 1.68), smallest ESS 5; 2,897 draws at the maximum tree depth. R-hat values near 2 mean the chains sit in different places: on 14 years of real NEE the main model (no cold limit) has more than one region of good fit, mostly in the leaf timing and leaf life parameters. The same setting passed on fake data (F23), so this is a property of the model and the real data, not of the code (Gelman et al. 2020, Section 5.1: a computational problem often points to a model problem).
- **Chain by chain (`diagnose_chains.py`, 2026-10-21): two separate modes.** Chains 2 and 3 agree (log posterior -4,228 and -4,227); chains 0 and 1 sit in clearly poorer places (48 and 69 lower; log-likelihood 37 and 46 worse). They differ mainly in the share of growth sent to the leaves (`f_fol` 0.15 to 0.17 against 0.065 to 0.07), LMA (about 270 against 330) and the starting soil pool (about 6,200 against 7,800). **No chain is believable:** Tier 2 plausible 0% in every chain, peak LAI about 1.0, `ceff` about 100, log prior -128 to -149 (only 1 prior draw in 1,000 is below -36.6). So the main model has both a model problem (implausible forest) and a sampling problem on top (separate modes, Section 5.6).
- Also submitted on 2026-10-21: the cold-limit model again with the D53 setting written out (Roihu job 2415129, `fit_model_error_cold_dense_ta95`, 72 h limit).

**F31. Cold-limit model, full fit with the D53 setting: fails D49 on divergences only; results not interpreted (Roihu job 2415129, dense mass matrix, target_accept 0.95, 13 h 35 min, finished 2026-10-22).**
- **7 divergences** of 4,000 draws (rule: none). Everything else passes clearly: largest R-hat 1.0055 (`ceff`), smallest ESS bulk 1,770 and tail 1,754, **no draw at the maximum tree depth** (against 3,969 of 4,000 with the diagonal setting, F29). The dense setting fixed the slow mixing; what is left is a few places where the step is still too large for the curvature.
- **Planned response, as written in D49 ("If it fails") and done in D53:** first look where the divergent draws sit, then raise `target_accept` (0.95 to 0.99), which makes the steps smaller. The rule of no divergences is kept as written.
- **Where the divergences sit (2026-10-23):** spread over all four chains (1, 2, 1 and 3), and the divergent draws have nearly the same values as the others (for example `ceff` 97.8 against 97.7, `cold_t_shutdown` -8.9 against -8.8); only a slight lean towards lower LMA and higher `c_lf` and `f_fol`. No chain and no extreme corner of one parameter stands out, so this is the case of steps slightly too large in a curved region, which a higher `target_accept` addresses; reparameterising (D27) is not indicated. These values were read only to locate the divergences, not as results.

**F32. The D63 ladder: fake-data checks and full fits (Roihu, all dense mass matrix with target_accept 0.99, finished by 2026-10-25). Only the D49 checks were read at this step.**

| Model | Free | Fake-data check | Truth inside 90% | Real fit (D49) | Fit time |
|---|---|---|---|---|---|
| M2 | 26 | PASS (F26, target_accept 0.95) | 24 of 26 | PASS: 0 divergences, R-hat 1.0026, ESS 1,949 | 27.3 h (job 2484088) |
| MA | 25 | PASS (job 2495588, 6.4 h) | 23 of 25 (need 21) | PASS: 0 divergences, R-hat 1.0026, ESS 1,756 | 22.0 h |
| MB | 24 | PASS (job 2495590, 7.4 h) | 22 of 24 (need 20) | PASS: 0 divergences, R-hat 1.0023, ESS 2,223 | 25.5 h |
| MC | 22 | PASS (job 2495592, 5.7 h) | 20 of 22 (need 18) | PASS: 0 divergences, R-hat 1.0016, ESS 1,708 | 27.0 h |
| M0 | 14 | **FAIL** (job 2495594, 3.9 h) | **9 of 14 (need 11)** | PASS: 0 divergences, R-hat 1.0038, ESS 1,084 | 4.9 h |

- **M2 now passes D49** with target_accept 0.99 (the 7 divergences of F31 are gone), as planned in F31. The M2 fit was submitted twice by mistake (jobs 2484088 and 2495596); both gave identical numbers, as expected from the fixed seeds.
- **MA, MB and MC** miss the same two quantities as F26 (`cr_onset` and `cold_t_full`). This is expected rather than reassuring: the free quantities of the true forest are the same draw each time, so these are not independent repeats.
- **M0 fails its fake-data check.** The sampler worked (all three computation checks pass), but the truth is outside the 90% interval for 5 of 14 quantities: `theta_roo`, `c_lf`, `ceff`, `c_fol_0`, `c_roo_0`. With 14 quantities at 90%, 5 or more misses happen by chance in about 1% of runs. The misses are canopy and root quantities, which suggests that, with everything else held, the posterior settles on a different canopy and root combination than the truth. Under the D63 rule (step 1), **M0 takes no part in the model choice**, although its real fit passes D49. It is reported with this failure.
- **M0 failure checked for a code error (2026-10-25): none found.** Every fixed quantity has exactly the same value in the truth and in the fitted model (for example `f_auto` 0.510134 and `f_fol` 0.0746952 in both). The sampler is healthy, and the M0 log density matched the hand computation in `tests/test_fix_sets.py`. The five misses are all near misses beside very narrow intervals (posterior width 8% to 28% of the prior's), and they go in compensating directions: more efficient leaves with less starting foliage (`ceff` 12.4 against a truth of 11.5; `c_fol_0` 453 against 491), and faster root turnover with a smaller starting root pool (`theta_roo` 0.0025 against 0.0016; `c_roo_0` 312 against 371). `c_lf` is missed by less than 0.001. So the posterior sits on a ridge just beside the truth, with intervals too narrow to reach it.
- **Two facts that help explain it; neither changes the verdict.** (1) M0's true forest is a different prior draw from the one used for MA to MC, because fixing 12 quantities changes which draw is the first plausible one, and it is an unusual draw (`cold_t_shutdown` -0.80 degC, above the top of its prior's 90% range, -2.03). (2) The true forest is drawn from the separate priors without the soft foliage rule (D44), while the fitted model includes that rule. With `f_fol` and `lma` fixed, M0 can only meet the rule by moving `c_fol_0`, `c_lf` and `ceff`, which are three of the five misses. This is a likely mechanism, not a tested one. Running M0 again with a different truth until it passes would be tuning the check to the answer, so it is not done. The D63 verdict stands.
- **Taking part in the choice (D63 step 2): M2, MA, MB and MC.**

**F33. Second check of the model code: the Edinburgh DALEC2 release (Williams 2019, doi:10.7488/ds/2504, `DALEC2 Feb19.f90`, read 2026-10-25). The model matches it; one small difference in day length, kept as is.**
- **What it is:** DALEC2 ("DALEC_CDEA") in Fortran 90, deposited by M. Williams in March 2019 as the code of López-Blanco et al. (2019). It is a second official copy, independent of the CARDAMOM C code (`DALEC_400`) used so far (D21 to D31).
- **Matches line by line:** all ten ACM coefficients and `psi_d`, `R_tot` (D22); conductance with half the temperature range; temperature on `Tmax` in ACM; respiration and decomposition on the midpoint temperature (D20); allocation in steps (foliage, then labile, then roots, wood the rest; same model as our Dirichlet form, D28); leaf fall anchored on the leaf-fall day `F_day`, not on its duration (D25, the typo in printed A8); leaf timing on the day of the run, not the day of year (D30); the offset polynomial in log(L - 1) and the lifespan 1.001 for the labile release (D31); wood litter straight to soil, root and leaf litter to the litter pool; LAI = C_fol / LMA; NEE = Ra + Rh - GPP. **So the model as its authors run it is confirmed by two separate official codes.** Re-read 2026-10-25 against the paper itself: of the "four corrections" carried over from the previous repository, only one changes an equation (A8 prints `cr_fall` where `d_fall` belongs; the paper's own sentence before A9 says the maximum leaf loss occurs at `t mod 365.25 = d_fall`, and both codes anchor on the leaf-fall day). The allocation one is an inconsistency inside Table 1 (the table entry for `f_woo` includes `f_roo`, the footnote below the table leaves it out). The two prose slips are not verified and do not affect the model; they are left out of the thesis. A9 does not matter: the paper says the authors use a sixth-order polynomial look-up for `psi`, which is the polynomial in both codes and in ours. In the thesis these are called differences between the printed equations and the official code, not errors in DALEC2.
- **One difference: the day-length formula.** This release uses declination `-asin(sin(23.45 deg) cos(2 pi (doy + 10) / 365))` and labels the formula in CARDAMOM `ACM.c` (declination `-23.4 deg cos(...)`, which is ours, D23) as "the old REFLEX project calculation but it is wrong". The two day-length formulas are otherwise the same (12 (1 + 2 asin(x) / pi) equals 24 acos(-x) / pi), so only the declination differs. **Size at FI-Hyy (61.85 N):** at most 4 minutes of day length (mean 2.3 minutes); the ACM day-length factor changes by -0.4% to +0.6% (spring and summer -0.4% to +0.1%; the largest relative changes are in late autumn and winter, when GPP is near zero).
- **Decision: kept as is; no refit.** The difference is far smaller than the measurement error and the fitted model error, and our formula is the one in the CARDAMOM code that DALEC_400 and the cold-limit code (DALEC_831) run with. Stated in the thesis as a known difference between official versions.

**F34. Test years and the model choice (D61, D63; Roihu job 2526482, 2026-10-25). By the rule fixed in D63, the chosen model is MA (the cold model with the starting wood pool fixed). No model gives a believable forest.**

| Model | Free | Test score (log predictive density, 810 days) | 90% coverage (all / DJF / MAM / JJA / SON) | RMSE | Bias (obs - model) | Tier 2 share | Peak LAI | Rh (g C m-2 yr-1) |
|---|---|---|---|---|---|---|---|---|
| M2 | 26 | -886.5 | 0.85 / 1.00 / 0.85 / 0.70 / 0.90 | 0.722 | -0.11 | 0% | 1.06 | 555 |
| **MA** | 25 | **-885.1** | 0.85 / 1.00 / 0.85 / 0.70 / 0.90 | 0.721 | -0.11 | 0% | 1.05 | 555 |
| MB | 24 | -887.3 | 0.85 / 1.00 / 0.86 / 0.69 / 0.91 | 0.723 | -0.12 | 0% | 1.05 | 555 |
| MC | 22 | -890.4 | 0.85 / 1.00 / 0.85 / 0.68 / 0.91 | 0.725 | -0.12 | 0% | 1.05 | 552 |
| M0 (not eligible, F32) | 14 | -992.1 | 0.83 / 1.00 / 0.89 / 0.61 / 0.89 | 0.794 | -0.21 | 5.8% | 1.45 | 530 |

- **D63 steps 2 to 4:** MA has the best score and is "clearly better" than each other eligible model by the D61 rule: by 1.5 over M2 (twice the widened standard error: 0.9), by 2.2 over MB (1.9) and by 5.3 over MC (4.3). No other model is within the margin, so step 4 (take the smallest of the equally good) has only MA to choose from.
- **Step 5 (safeguard) does not apply:** no eligible model is believable (Tier 2 share 0% for all four), so there is no believable alternative to prefer. The rule picks **MA**.
- **Practical size, stated honestly:** the differences are statistically clear but practically tiny. 1.5 log units over 810 days is 0.002 per day; RMSE differs by under 0.5% and coverage by at most 2 points between M2, MA, MB and MC. The standard errors are small because the models are nearly the same and make nearly the same errors on the same days. The widening uses only the lag-1 correlation of consecutive days; if the errors are correlated over longer spans (for example by season), the standard errors are still too small and "clearly better" is less certain than the rule says. Reading: **fixing the starting wood pool costs nothing and may help a little; fixing more starts to cost (MB, MC slightly, M0 a lot: 107 log units worse).** The scores fall steadily as more is fixed: -885.1, -887.3, -890.4, -992.1.
- **PSIS-LOO on the fitted years (reported only):** all four within 2 of each other (-3,746 to -3,748, standard error about 77; differences 0.2 to 1.4 with standard errors 0.3 to 1.4), so the fitted years cannot tell them apart either. Pareto k fine for every day.
- **The forest is still not believable in any model:** peak LAI about 1.05 (site 2.5 to 4), Rh about 555 g C m-2 yr-1 (site 290 to 370), Tier 2 share 0%. The cold limit fixed the spring timing but not the thin canopy. This matches F22 and F28: with the ACM light-use efficiency `a7` fixed, a thin canopy is the model's only way to lower its light-use efficiency. This is a main result for RQ1 and RQ3: NEE alone gives a model that predicts NEE reasonably well (coverage 85%) from a forest whose inside is wrong.
- **Shared weaknesses of all models on the test years:** coverage 100% in winter and only 68% to 70% in summer (one model-error size for the whole year is too wide in winter and too narrow in summer); bias about -0.11 g C m-2 d-1 (the model releases a little more carbon than measured), so the yearly sink is underestimated (2011: about -75 against -127 observed; 2013: about -272 against -313).
- **Found while reading:** all 810 test days are in 2011 to 2013. 2014 has no day that passes the D12 rule (no sums are printed for 2014), so D61's "2011 to 2014" is in practice 2011 to 2013. To check in the data and state in the thesis.
- Minor: the evaluation csv and figure file names list only the first three fits (the code shortens names of more than three); the files contain all four.

**F35. The literature values (D64): the only independent DALEC2 posterior for this location describes a different ecosystem from the stand (CARDAMOM 2001-2010, Bloom et al. 2016; cell centred at 61.5 N, 24.5 E; read 2026-10-27 with `scripts/07_literature_priors.py`).**

| Quantity | CARDAMOM cell (median) | Site prior median (D34 to D39) | Site measurement |
|---|---|---|---|
| LMA | **22.8 g C m-2** | 185 | about 165 to 210 (D38) |
| foliage pool (2001-2010 mean) | **28 g C m-2** | 590 (`c_fol_0`) | about 590 (D38) |
| foliage residence time | **0.40 yr** | needle life 3.8 yr (`c_lf` 0.262) | 3 to 5 yr |
| mean LAI (foliage / LMA) | 1.2 | about 3.2 | 2.5 to 4 (peak) |
| share of GPP to foliage | 0.011 | 0.075 | |
| `f_auto` | 0.50 | 0.51 | |
| wood pool | 10,760 g C m-2 | 5,900 | 5,900 to 6,300 |
| soil pool | **23,610 g C m-2** | 6,100 | 6,560 |
| wood residence time | 30.7 yr | 51 yr | |
| soil residence time (at the site's temperatures) | 54.1 yr | 32 yr | |
| GPP | 1,068 g C m-2 yr-1 | | 952 to 1,104 |
| Rh | 538 g C m-2 yr-1 | | 290 to 370 |

- **Not a quirk of one cell:** across southern and central Finland (60 to 64 N, 22 to 30 E) the CARDAMOM cells have LMA about 23 (median of cells) and a foliage residence time of 0.35 to 0.42 yr, and soil carbon 20,000 to 50,000 g C m-2. CARDAMOM represents this region with a thin canopy that is renewed every year (deciduous-like), fitted to 1 degree MODIS LAI, and with the soil carbon of a global soil map (which includes peat).
- **What it agrees with, and why this matters (RQ3, Discussion):** with the same fixed ACM coefficients (D22), CARDAMOM reproduces the site's GPP (1,068) with a mean LAI of about 1.2, and its Rh (538) is above the site range, much like our NEE fit (peak LAI about 1.05, Rh about 555, F34). Two independent DALEC2 analyses, one constrained by satellite LAI and soil carbon, the other by tower NEE, both end with a thin canopy and high Rh at this location. This supports the reading of F22 and F28 (a property of DALEC2 with its fixed ACM coefficients at a boreal site) over an artefact of our likelihood. It is support, not proof: a 1 degree cell is not the stand, MODIS LAI includes lakes and other land, and a mean LAI is not a peak LAI.
- **Why these values do not work as priors for the stand, by the D64 rules as written:**
  1. **Their centres conflict with the site measurements** by factors of 4 to 21 (LMA 8 times lower, foliage pool 21 times smaller, needle life about 10 times shorter, soil carbon 4 times more).
  2. **After doubling, the spreads are very wide** (log sd 1.5 to 2.3 for wood and soil turnover, LMA and the wood pool: a factor of 4.5 to 10 per sd), so they carry little information, which is the opposite of what the supervisor's idea needs.
  3. **The four Dirichlet weights are all below 1** (0.029, 0.006, 0.022, 0.201; sum 0.26). Below 1, the density becomes infinite as a share approaches 0, which breaks the rule that every prior is smooth with no walls (D33) and is hard for NUTS.
  4. **A dry run of the prior predictive check (made-up weather, 2026-10-27, code check only, not reportable):** 85% of draws physically impossible (LAI above 12 in 85%), 0.7% plausible. Combining CARDAMOM's LMA (23) with the site's foliage prior gives impossible canopies. The real check on the site weather is still to run (D65).

**F36. The D63 ladder with the D66/D67 code (direct `theta_min` prior over the published range, allocation in beta steps; Roihu jobs 2684863 to 2684873, finished 2026-11-01). Only the pass/fail lines were read.**
- **Prior predictive check of the cold model: PASS** (Tier 1 3.5%, Tier 2 57.7%, with the soft foliage rule).

| Model | Free | Fake-data check: computation | Truth inside 90% | Real fit (D49) |
|---|---|---|---|---|
| M2 | 25 | pass (largest R-hat 1.0098, `theta_min`; ESS 482) | **20 of 25, need 21: FAIL** (f_auto, frac_roo, cr_onset, lma, ceff) | pass |
| MA | 24 | **FAIL** (R-hat 1.0112, `theta_min`; ESS 344) | 20 of 24, need 20 | pass |
| MB | 23 | pass | **18 of 23, need 19: FAIL** (same five as M2) | **FAIL: 1 divergence** |
| MC | 21 | pass | 18 of 21, need 17: pass | **FAIL: 1 divergence** |
| M0 | 14 | pass | **8 of 14, need 11: FAIL** | pass |

- **Under the D63 rule, no model is eligible** (none passes both its fake-data check and its fit).
- **Compared with F32 (old prior via `h`, Dirichlet form), where M2, MA, MB and MC all passed:** two things changed at once. (1) `theta_min` now has a prior about 10 times wider, and it is the slowest-mixing quantity in four of the ten runs (largest R-hat in the M2, MA and MC fake-data checks and the MB fit), and the smallest effective sample sizes of the fake-data checks for M2 to MC fell from about 1,100 to 2,300 (F26, F32) to 344 to 781. (2) The true forest of the fake-data checks is a different prior draw, because the priors changed and the random draws with them; the new misses are canopy and allocation quantities (`lma`, `ceff`, `f_auto`, `frac_roo`, `cr_onset`), not `theta_min`. The two causes are not yet separated: next step is to look at where the truth sits and how `theta_min` trades off with other quantities, before anything is changed (D49, "if it fails").

**F37. The D68 round: two true forests, 2,000 draws per chain, target_accept 0.995 (Roihu jobs 2709839 to 2709850, finished 2026-11-03/04). Under D63, MC is the only eligible model; the test-year scores (job 2804864) are still to be read.**

| Model | Free | Fake-data computation (both truths) | Truth inside, truth 1 / truth 2 | Pooled (rule) | Real fit (D49) | Inside 90% band, 1997-2010 |
|---|---|---|---|---|---|---|
| M2 | 25 | pass (R-hat at most 1.0069, ESS at least 851) | 20 / 20 of 25 | **40 of 50 (need 42): FAIL** | pass (0 div., R-hat 1.0014, ESS 2,379) | 89.1% |
| MA | 24 | pass (1.0049, 773) | 19 / 20 of 24 | **39 of 48 (need 41): FAIL** | pass (0, 1.0013, 3,271) | 89.0% |
| MB | 23 | pass (1.0055, 1,010) | 18 / 18 of 23 | **36 of 46 (need 39): FAIL** | pass (0, 1.0017, 2,845) | 89.3% |
| MC | 21 | pass (1.0023, 833) | 17 / 18 of 21 | **35 of 42 (need 35): PASS** | pass (0, 1.0031, 2,158) | 89.1% |

- **Computation is now clean everywhere:** no divergences in any of the 12 runs, all R-hat below 1.007, all ESS above 770. The D68 sampler settings solved the F36 computation failures. Every fake-data failure is on the truth criterion only.
- **The misses come in the same groups in every model**, not at random:
  - truth 1: `f_auto`, `lma`, `ceff`, `cr_onset` (all four models), plus `frac_roo` (M2, MA, MB);
  - truth 2: `cr_onset`, `ceff`, `cold_t_full` (all four models), plus `theta_lit` (M2, MA, MB), `c_lit_0` (M2) and `theta_woo` (MB).
  `ceff` and `cr_onset` are missed for both truths in every model.
- **Reading (to be checked with the correlations, not used to change any rule):** the groups match trade-offs expected from the equations. `f_auto`, `lma` and `ceff` act together on photosynthesis seen by NEE (leaf area = foliage / LMA; NEE sees (1 - f_auto) x GPP); `cr_onset` and `cold_t_full` both set the start of spring. When quantities trade off, they miss together. The D63 misses rule counts misses as if they were independent (binomial), so for correlated quantities it is stricter than intended: a pass or fail can turn on one shared trade-off. This is reported as a limitation of the check, not used to change the result (D68: the rules are not changed again in this round). A full calibration check would need many true forests (simulation-based calibration), which is out of reach here.
- **Why MC passes and MB does not:** MC also fixes the soil and wood turnover rates, and with them the `theta_lit` (and in MB `theta_woo`) misses disappear, consistent with a trade-off along the litter-to-soil path. MC passes with no margin (35 of 35 needed).
- **Real fits:** the four models fit 1997-2010 almost identically (89.0 to 89.3% inside the 90% band). In every model the band is too wide in winter (DJF 99.9%) and too narrow in summer (JJA 75.4 to 75.9%), with spring 89.7 to 89.9% and autumn 92.2 to 92.5%; average misses are small (-0.09 to +0.03 g C m-2 d-1). One model error `sigma_m` for the whole year is too small for summer and too large for winter: a limitation of the error model (L1), reported, not changed.
- **Consequence under D63:** only MC passes both its fake-data check and its fit, so MC is the only model that takes part in the choice. The test-year scores of all four are still computed and reported (D63 point 6), with MC's believability (D40 Tier 2).

**F38. Test years and the D63 choice for the D68 round (Roihu job 2804864, 2026-11-04). By the D63 rule the chosen model is MC (the cold model with the starting wood and soil pools and the wood and soil turnover rates fixed). The four models predict almost the same; none gives a believable forest.**

| Model | Log predictive density, 810 test days | Behind M2 (D61 rule) | Coverage 90%: all / DJF / MAM / JJA / SON | RMSE | Bias (obs. - model) | Tier 2 share | Peak LAI | GPP | Rh |
|---|---|---|---|---|---|---|---|---|---|
| M2 | -885.6 | best | 85.2 / 100 / 85.7 / 69.6 / 90.0% | 0.7217 | -0.111 | 0 | 1.05 | 1,076 | 556 |
| MA | -886.4 | 0.7 (limit 0.6): clearly behind | 84.8 / 100 / 85.2 / 68.8 / 90.0% | 0.7219 | -0.111 | 0 | 1.05 | 1,079 | 555 |
| MB | -886.6 | 1.0 (limit 0.8): clearly behind | 84.6 / 100 / 84.6 / 68.4 / 90.0% | 0.7222 | -0.112 | 0 | 1.05 | 1,077 | 555 |
| MC | -888.6 | 3.0 (limit 2.8): clearly behind | 84.9 / 100 / 85.2 / 68.8 / 90.4% | 0.7234 | -0.115 | 0 | 1.05 | 1,076 | 553 |

(Units: NEE g C m-2 d-1; GPP and Rh g C m-2 yr-1, medians over 1997-2010; site peak LAI 2.5 to 4, site Rh 290 to 370.) Test-year sums, observed / model median (all four models within 1 g C m-2): 2011 -127 / -76; 2012 -212 / -215; 2013 -313 / -272. PSIS-LOO on 1997-2010 (reported only): all four -3,747 to -3,748, differences within 1 standard error.

- **The D63 choice is MC.** Only MC passes both its fake-data check and its fit (F37), so it is the only model that takes part; points 3 to 5 of D63 have nothing to compare. The safeguard (point 5) does not apply: no model is believable (Tier 2 share 0 for all four).
- **M2 scores "clearly better" than every reduced model under the D61 rule, but by very little.** The largest difference (M2 against MC) is 3.0 log units over 810 days, about 0.004 per day; RMSE differs in the third decimal (0.7217 against 0.7234) and the annual sums agree within 1 g C m-2. The standard errors are small because the models are nearly the same and are compared day by day, so differences too small to matter in practice still pass the 2-standard-error rule. This is reported as a limitation of the rule (statistical, not practical, difference), not used to change the choice. The model that predicts best (M2) is not eligible because it failed its fake-data check.
- **What fixing cost:** holding the two slow pools and their turnover rates fixed (21 instead of 25 quantities) changed the fit to 1997-2010 and the predictions for 2011-2013 by amounts that are negligible in practice. For NEE at this site, these four quantities are close to irrelevant: the simpler model represents the NEE as well as the full one.
- **The forest is not believable in any model, and in the same way:** peak LAI about 1.05 (site 2.5 to 4) and Rh about 555 (site 290 to 370) in all four. Since fixing wood and soil changes neither, the cause is not in those quantities; it points to the canopy (foliage and LMA) and to how respiration is split. This is the main item for MC's posterior analysis and for D65.
- **Test years:** coverage of the 90% band 85% overall, 100% in winter and about 69% in summer, the same pattern as the fitted years (F37) but stronger. The model takes up less carbon than observed on the test days (bias -0.11), mostly in 2011 and 2013.

**F39. The chosen model's posterior (MC, `fit_model_error_MC_cold_dense_ta995_d2000`; summary table from `talk_figures.py`, read 2026-11-04). The same thin, efficient canopy as F21 and F22: 9 of the 21 quantities have their posterior mean outside the prior's 90% range.**

| Quantity | Prior 5% to 95% | Posterior mean (90%) | Where the mean sits in the prior |
|---|---|---|---|
| `ceff` | 5.6 to 17.9 | 98 (89 to 108) | far above (quantile 1.0000) |
| `lma` (g C m-2) | 123 to 279 | 442 (359 to 532) | above (0.9998) |
| `c_fol_0` (g C m-2) | 360 to 966 | 159 (120 to 202) | far below (0.0000) |
| `c_roo_0` (g C m-2) | 145 to 390 | 670 (423 to 958) | above (0.9997) |
| `c_lf` (yr-1) | 0.21 to 0.32 | 0.37 (0.35 to 0.40) | above (0.9989): needles live about 2.7 years |
| `temperature_exponent` | 0.019 to 0.060 | 0.075 (0.073 to 0.077) | above (0.989): Q10 about 2.1 |
| `f_auto` | 0.39 to 0.63 | 0.30 (0.25 to 0.35) | below (0.0009) |
| `theta_lit` (d-1) | 6.2e-4 to 2.4e-3 | 4.8e-4 (3.6e-4 to 6.2e-4) | below (0.012) |
| `c_lab_0` (g C m-2) | 34 to 175 | 189 (145 to 239) | above (0.964) |

- Inside the prior's 90% range: the three step fractions, `theta_roo`, `theta_min` (low end, quantile 0.08), the four phenology dates and rates, `c_lit_0`, and the two cold thresholds (`cold_t_full` 1.6 degC, low end of its range). `sigma_m` 0.52 g C m-2 d-1.
- **Reading:** the posterior is the thin canopy of very efficient leaves found in F21 and F22 (peak LAI about 1.05, `ceff` about 100, LMA about 440), with less plant respiration (`f_auto` 0.30) and more decay (Rh about 553, site 290 to 370) and a steep temperature response. F22 showed that NEE strongly prefers this canopy because, with `a7` fixed at the Loobos value, it is the model's only way to make GPP follow sunlight; F28 showed that estimating `a7` moves the canopy towards believable values. Fixing the slow pools (MC) did not change any of this (F38).
- **Computation:** largest R-hat 1.0031 (`frac_lab`), smallest ESS 2,158; every quantity passes D49.
- **Next:** the prior sensitivity (D69, job submitted 2026-11-04) is expected to flag these quantities as prior-data conflicts; it is read as stated in D69, not used to change any prior.

**F40. Prior sensitivity of MC, all priors scaled together (D69; Roihu job 2805470, 2026-11-04; log prior and log likelihood computed in 0.9 min).** Threshold 0.05 (ArviZ default).

| Reading | Quantities (prior sensitivity / likelihood sensitivity) |
|---|---|
| Data decide, prior hardly matters | `theta_roo` (0.039 / 0.142), `temperature_exponent` (0.031 / 0.085), `ceff` (0.035 / 0.093) |
| Prior decides, data weak | `frac_roo` (0.083 / 0.022) |
| Prior and data pull against each other ("potential prior-data conflict") | the other 18, strongest: `c_roo_0` (0.605 / 0.529), `theta_lit` (0.526 / 0.570), `c_lab_0` (0.306 / 0.202), `c_fol_0` (0.278 / 0.411), `c_lf` (0.251 / 0.205), `c_lit_0` (0.243 / 0.150), `d_onset` (0.233 / 0.292), `frac_lab` (0.220 / 0.083) |

- **Two kinds of disagreement with the prior.** `ceff` and `temperature_exponent` sit far outside their priors (F39) yet have low prior sensitivity in this joint table. **Corrected by F41 for `ceff`:** scaled on its own, its prior does matter (0.192); the low joint value came from all priors moving together. For `temperature_exponent` F41 confirms that the data decide. Quantities such as `c_roo_0`, `theta_lit`, `c_fol_0` and `c_lf` are also outside or at the edge of their priors but with high prior sensitivity: there the prior is still holding them back, and the answer depends on it. (Power-scaling measures how much a small change of the prior moves the posterior, not how far the posterior is from the prior.)
- **Caution:** the joint scaling flags 18 of 22 quantities, including some well inside their priors (`d_onset`, `d_fall`, `cr_onset`). Because all priors are scaled at once, a quantity correlated with a conflicting one (for example through the canopy and pool trade-offs) moves with it. The own-prior table (D69 addition) is run to separate these.

**F41. Prior sensitivity of MC, each quantity's own prior scaled on its own (D69 addition; 2026-11-04, from the saved log prior and log likelihood).** Threshold 0.05.

| Reading | Quantities (own prior / likelihood) |
|---|---|
| **Data decide; own prior hardly matters** (9) | `temperature_exponent` (0.015 / 0.085), `d_onset` (0.030 / 0.292), `cr_onset` (0.026 / 0.252), `d_fall` (0.015 / 0.144), both cold thresholds (0.004 / 0.072; 0.023 / 0.166), `theta_min` (0.044 / 0.193), `c_lit_0` (0.043 / 0.150), `sigma_m` (0.000 / 0.159) |
| **Own prior decides, data weak** (1) | `frac_roo` (0.066 / 0.022) |
| **Own prior and data pull against each other** (12) | strongest: `c_roo_0` (0.533 / 0.529), `c_fol_0` (0.343 / 0.411), `lma` (0.283 / 0.061), `f_auto` (0.265 / 0.156), `c_lf` (0.213 / 0.205), `ceff` (0.192 / 0.093), `theta_lit` (0.125 / 0.570), `c_lab_0` (0.118 / 0.202), `frac_fol` (0.114 / 0.158); near the threshold: `cr_fall` (0.087), `frac_lab` (0.078), `theta_roo` (0.065) |

- **The joint flags on the phenology dates were a side effect** of scaling all priors together (F40): on their own, `d_onset`, `cr_onset` and `d_fall` are decided by the data.
- **Clear split for RQ1.** NEE decides the timing of the leaf season, the temperature response of respiration, the cold thresholds, the litter-to-soil transfer (`theta_min`, low end of its range), the starting litter pool and the model error. In the canopy and allocation group (`ceff`, `lma`, `c_lf`, `f_auto`, `frac_fol`, starting foliage, labile and root pools, litter decay) the posterior depends on the prior as well as on the data, and the data push most of them outside their priors (F39). These are the quantities of the thin-canopy solution (F22): there the posterior is a compromise between the prior and a model that cannot make GPP follow light with a normal canopy, not a clean estimate of the forest.
- **Correction to F40:** `ceff` is not "decided by the data"; its own prior matters.
- Power-scaling is local (a small change of the prior) and uses re-weighted draws; it describes the fitted posterior and changes no prior and no choice (D69).

**F42. RQ1 summary for MC: width ratio, position in the prior and own-prior sensitivity together, and fake against real width ratios (`scripts/12_rq1_summary.py`; fake data: mean of truths 1 and 2, 2026-11-04).** Table: `figures/rq1_fit_model_error_MC_cold_dense_ta995_d2000.csv`; figure AA.

- **Three groups (own-prior reading, F41):**
  - *NEE decides* (9): the temperature response (width ratio 0.04), both cold thresholds (0.10, 0.16), the start and length of leaf growth (`d_onset` 0.20, `cr_onset` 0.16), leaf fall (`d_fall` 0.31), `theta_min` (0.56), `c_lit_0` (0.52) and `sigma_m`. Only the temperature exponent and `cold_t_full` sit outside their priors (position 1.40 and -1.15).
  - *Prior and NEE pull against each other* (11, plus `frac_lab` near the threshold): the canopy and allocation group. Most sit outside their priors: `ceff` (position 3.9), `c_fol_0` (-2.7), `lma` (2.1), `c_roo_0` (2.1), `f_auto` (-1.9), `c_lf` (1.8), `theta_lit` (-1.4), `c_lab_0` (1.1).
  - *The prior decides* (1): `frac_roo`.
- **Fake against real (mean of the two true forests; Roihu run 2026-11-04).** 12 of 21 quantities are within 0.1 of the diagonal: on real data NEE taught about as much as the fake-data checks said it could. Real data taught clearly less for `theta_min` (0.56 against 0.19), `c_lit_0` (0.52 against 0.23), `c_roo_0` (0.83 against 0.62), `cold_t_full` (0.16 against 0.02), `c_fol_0` (0.53 against 0.40), `theta_lit` (0.40 against 0.28) and `f_auto` (0.50 against 0.40); part of this is expected, because the real fits also carry the model error `sigma_m` and the fake ones do not. Real data taught more for `cr_onset` (0.16 against 0.43) and `c_lf` (0.39 against 0.49), as with one truth; for `c_lf`, which sits outside its prior, a narrower posterior than the design allows suggests it is pinned by the model-data mismatch rather than learned.

**F43. RQ2, where MC's predictive uncertainty comes from (figure X, `talk_figures.py --predict`, 200 posterior draws; 2026-11-04).** Standard deviations by month, g C m-2 d-1:
- **Model error `sigma_m`: 0.52 in every month**, the largest source all year.
- **Measurement error (mean RANDUNC on the days used):** about 0.06 in winter, rising to 0.27 in July.
- **Parameter uncertainty (spread of daily NEE across posterior draws): 0.01 to 0.05**, largest in May, the smallest source in every month.
- **Reading.** The predictive band is set almost entirely by the model error. Because `sigma_m` is one value for the year, the total band hardly changes with the season (about 0.52 in winter, about 0.59 in July: the square root of the sum of the squares), while the real misfit is much larger in summer: hence the 70 to 76% summer and 100% winter coverage (F37, F38). Parameter uncertainty is small even where the parameters disagree with their priors (F39 to F42): the posterior is confident about a forest that is not believable. More NEE data would narrow the smallest source; a better model structure (or a model error that changes with the season) would address the largest.

**F44. RQ3, which quantities trade off in MC (figure U, the six most correlated pairs on the prior's scale; 2026-11-04).**

| Pair | Correlation | Why, from the model equations |
|---|---|---|
| LMA and `c_lab_0` | +0.77 | more starting sugar store means more foliage in the spring flush; a heavier leaf (higher LMA) keeps leaf area, which NEE sees, the same |
| `theta_lit` and `c_roo_0` | -0.75 | dying roots feed the litter pool; more root carbon at the start needs slower litter decay for the same decay respiration |
| `frac_lab` (s_lab) and `c_lab_0` | +0.73 | both set how much carbon passes through the labile store to the leaves |
| LMA and `c_fol_0` | +0.71 | leaf area = foliage carbon / LMA: NEE fixes their ratio, not each |
| `frac_lab` (s_lab) and LMA | +0.71 | more carbon to the labile store, more foliage, so a heavier leaf for the same leaf area |
| `f_auto` and `ceff` | +0.65 | NEE sees (1 - `f_auto`) x GPP: more plant respiration needs more photosynthesis (predicted before the results were read, 2026-11-03) |

- **Four of the six pairs are one trade-off:** {LMA, `c_fol_0`, `c_lab_0`, `frac_lab`}. NEE fixes how much leaf area the forest has, but not how much carbon is in the leaves or how heavy each leaf is. This is the "thin canopy" group (F22, F39).
- **Every quantity in these pairs is in the "prior and data pull against each other" group of F41.** Where NEE fixes only a combination, the prior decides where along it the posterior sits; that is why their priors still matter. RQ1 and RQ3 explain each other.
- **The fake-data misses group the same way** (F37: `f_auto`, `lma`, `ceff` missed together).
- Expected pairs not among the six: `theta_som` and `c_som_0` (both fixed in MC, so no correlation exists); `theta_min` and `theta_lit`; `ceff` and LMA; `cr_onset` and `cold_t_full`. Their values are read from the full correlation matrix (all pairs with |r| of 0.5 or more), still to be listed.
- **What would separate them (for the Discussion):** leaf area or needle-mass measurements (LMA, foliage), root biomass (roots and litter), and a measured split of respiration (autotrophic against heterotrophic, `f_auto`).

## Part 5. Open items

### Still open

**OPEN-3. Spring onset number.** The previous EDA found a median observed onset of day 96 (range 58 to 117) with the 10-day rule, while project notes mention an observed onset of day 110 in the model comparison. The new code gives a median of day 100 (range 66 to 137) for 1997-2010 (F10), so three numbers now exist. These probably come from different definitions or data. Pick one definition (D16 uses the 10-day rule) and use it consistently in the thesis.

**OPEN-6. Small checks still to do** (none blocks the next step):
- D31: read the printed equation A9 in the Bloom and Williams (2015) PDF. **Narrowed 2026-10-16:** a text extraction of the PDF reads A9 as `2 sqrt(pi) log(1 - c_lf) psi - exp(-psi^2) = 0`, i.e. the previous repository's form. Text extraction loses fraction bars, so one look at the typeset PDF settles it: if `log(1 - c_lf)` sits *under* `2 sqrt(pi) psi` as a fraction, the paper matches the official code and the previous repository mis-copied it; if it is a product on the same line, the paper has a typo. Either way the model follows the code (D31).
- D37: confirm the leaf onset timing (Schiestl-Aalto et al. 2013 or Schiestl-Aalto 2017) and the needle fall timing (monthly litterfall data).
- F9: confirm the length of the frozen-soil season at Hyytiala from the soil temperature data (the previous notes say about four to six months).
- ~~D19: correct `research_plan_v2.md`, which still lists `VPD_F` as a driver.~~ Done 2026-10-16: the driver paragraph now lists the D19 inputs and says why VPD and precipitation are not used.
- D2: add `pymc`, `pytest` and `netcdf4` to `environment.yml`; freeze `environment.lock.yml` before producing thesis results.
- D44: raise with the supervisor that "no EDCs" became "a soft EDC on the foliage pool".
- ~~D55: check the official cold-limit code and its priors.~~ Done 2026-10-16: found in DALEC_831 (CARDAMOM v2.2); see D55.
- References marked "check" below: verify against the original papers.

**OPEN-7. Choice of sampler.** NUTS as the only sampler is carried over from the previous repository. Its settings are now fixed (D49); a dense mass matrix is being tested (D52); other alternatives are not yet re-assessed. Decide after the fake-data run: if it is too slow on the author's machine, the options are running on Roihu, or nutpie (another NUTS implementation that PyMC can call, often faster and with better warm-up, `nuts_sampler="nutpie"`).

**OPEN-8. How to respond to F17 (real NEE pulls the parameters out of the believable region).** Options, all already recorded as planned changes, none tuned to the result: (a) add model error to sigma (L1 option C, Gelman et al. 2020, Sections 7.1 and 10.5), so the likelihood stops treating DALEC2 as able to match each day within RANDUNC; (b) the frost switch (D24), aimed at the spring misfit; (c) both; (d) a Laplace error (L2). Plan: quick fits of (a), (b) and (c) first, each judged with `diagnose_chains.py` (are the chains believable and do they agree?). **Status 2026-10-16 (F18):** (a) removes the separate modes and is now the main error model (L1, changed 2026-10-16); (b) alone does nothing. The remaining implausible canopy leaves three routes (Gelman et al. 2020, Sections 5.2, 6.1, 7): tolerate and report it; simplify (fix the canopy parameters at site values); or expand (a cold limit on photosynthesis, Stettz et al. 2022). To discuss with the supervisor. Gelman et al. (2020, Section 5.9) warn against "haphazardly tweaking the model until computational problems disappear", so every option must have its own reason, and all results are reported, not only the one that works. **Status 2026-10-23:** the main model is not fitted again and the cold-limit model is the base (D62); the simplify route is now the ladder of reduced cold models (D63).

**OPEN-9. Literature posterior values as fixed values or priors (supervisor's suggestion, 2026-10-26).** The idea: take fitted (posterior) DALEC parameter values from similar studies and use them as fixed values or priors, starting from the simplest model. **Literature search, 2026-10-26:**

| Source | Model | Site / data | Numbers published? | Data separate from ours? | Use |
|---|---|---|---|---|---|
| Bloom et al. (2016), CARDAMOM 2001-2010 global, doi:10.7488/ds/316 (CC-BY 4.0) | DALEC2 (`DALEC_400`, our code) | 1 degree grid; MODIS LAI, HWSD soil carbon; **no NEE** | yes, per grid cell: median, spread and percentiles of allocation shares (of GPP), residence times of foliage, roots, wood, litter, soil, the six pools, LMA | **yes** | **best candidate for priors** (FI-Hyy grid cell); units to confirm (no README) |
| Exbrayat and Williams (2018), pan-Arctic CARDAMOM 2000-2015, doi:10.7488/ds/2334; López-Blanco et al. (2019) | DALEC2 (Williams 2019 code) | 1 degree; MODIS LAI, NCSCD soil carbon, biomass; no NEE | allocation, pools, transit times (median and intervals); parameter posteriors only as a figure (Fig. S2) | yes | possible, if the domain covers 61.8 N (not confirmed; the only Finnish tower used is Sodankylä) |
| Stettz et al. (2022) | DALEC2cold (CARDAMOM) | Niwot Ridge, US; tower GPP 2000-2018 | cold thresholds only: T0 = -7.1 +/- 1.1 degC, Tg = 6.0 +/- 2.6 degC | yes (other site) | prior for the two thresholds (widened for transfer) |
| Famiglietti et al. (2021) | 16 DALEC variants (C1: 23 parameters, CDEA phenology, ACM1, like DALEC2; S1: 11, E1: 17) | **includes FI-Hyy**; first 5 years (about 1999-2003) of monthly NEE, LAI, biomass | **no** parameter values (only skill metrics, figshare) | **no** (FI-Hyy NEE overlaps ours) | comparison only; could ask the authors; their model-complexity ladder is the closest published analogue to D63 |
| Fox et al. (2009), REFLEX | original DALEC (5 pools), not DALEC2 | Loobos (NL) and Hesse (FR) | no table of posteriors; priors only (Table 3) | yes | already used for the `ceff` prior range (D42) |
| Bloom and Williams (2015) | DALEC2 | three temperate AmeriFlux sites; LAI and soil carbon, NEE only to validate | priors in Table 1 | yes | already used |

**Rules for using them (to agree with the supervisor before any run):** (1) the other study's data must be separate from ours, otherwise the same NEE is counted twice (rules out Famiglietti et al. 2021 at FI-Hyy as a prior); (2) same model version and the same meaning of each parameter (for example residence time against rate at 0 degC; share of GPP against share of NPP; LAI convention); (3) a value from another site or a 1 degree cell is widened before use; (4) fixing at a literature value removes its uncertainty and carries over any bias, so an informative prior is preferred to fixing. **Already true:** our priors are based on FI-Hyy measurements (D34 to D39), and D63 fixed quantities at those site-based medians. **Already tested:** fixing canopy quantities at literature values made the fit much worse (F20, F22: `ceff` at the REFLEX median 10 gave LAI 8 to 14 and a log-likelihood about 110 to 170 lower), and M0 (12 quantities fixed) failed its fake-data check (F32) and scored far worse (F34). Status 2026-10-27: discussed; the plan is D64 (literature priors as priors, a ladder from the simplest model).

### Planned sensitivity runs (decided in principle, to run after the main fit)

| Run | What changes | Entry |
|---|---|---|
| Frost switch | GPP stopped when the midpoint temperature is below -2 C | D24 |
| RANDUNC only (the original main model) | sigma_t = RANDUNC_t, no model error | L1 (main model changed 2026-10-16) |
| Wider `f_auto` prior | sd 0.10 instead of 0.07, only if the posterior is about as wide as the prior | D34 |
| Wider priors where NEE gave no information | only for parameters whose posterior is about as wide as the prior | D35, way forward 3 |
| Day-to-day correlated errors | for example AR(1) errors, if time allows | F11 |

### Closed

**OPEN-1. Evaluation mask.** Closed 2026-10-21 (D61): the test uses the fitting rule (810 days); RMSE and bias are also reported on the 1,418 QC-only days.

**OPEN-2. Which days each sigma figure uses.** Closed 2026-10-08 (D15).

**OPEN-4. Daily temperature extremes (F7).** Closed 2026-10-12: derived from the half-hourly file (D18), confirmed on the real files (F10).

**OPEN-5. Model structure.** Closed 2026-10-13: M1 to M8 and M10 decided (D21 to D28, D30); M9 moved to the priors (D29) and closed there by D38.

### Carried over from the previous repository

Re-assessed (2026-10-12 and 2026-10-13):
- Allocation (D28), ACM form (D21), LAI convention (D29, D38), `d_fall` (D25): changed or confirmed.
- No VPD driver: kept, with the literature check in D19.
- No EDC constraints: hard accept/reject rules are still excluded (they break NUTS), but a smooth steady-state version on the foliage pool was added (D44).

Not yet re-assessed:
- NUTS as the only sampler (OPEN-7).

---

## References

- Bacour, C. et al. (2023). Assimilation of multiple datasets results in large differences in regional- to global-scale NEE and GPP budgets simulated by a terrestrial biosphere model. Biogeosciences 20, 1089-1111.
- Bloom, A. A. and Williams, M. (2015). Constraining ecosystem carbon dynamics in a data-limited world: integrating ecological "common sense" in a model-data fusion framework. Biogeosciences 12, 1299-1315.
- Bloom, A. A., Exbrayat, J.-F., van der Velde, I. R., Feng, L. and Williams, M. (2016). The decadal state of the terrestrial carbon cycle: global retrievals of terrestrial carbon allocation, pools, and residence times. PNAS 113(5), 1285-1290. (Version of DALEC2 in the CARDAMOM code `DALEC_400`.)
- Brunner, I. et al. (2013). Fine-root turnover rates of European forests revisited: an analysis of data from sequential coring and ingrowth cores. Plant and Soil 362, 357-372.
- CARDAMOM framework, source code. github.com/CARDAMOM-framework/CARDAMOM, files `C/projects/CARDAMOM_MODELS/DALEC/DALEC_400/DALEC_400.c`, `PARS_INFO_400.c`, `DALEC_ALL/ACM.c`, `DALEC_ALL/offset.c` (read 2026-10-13).
- Chuter, A. M., Aston, P. J., Skeldon, A. C. and Roulstone, I. (2015). A dynamical systems analysis of the data assimilation linked ecosystem carbon (DALEC) models. Chaos 25(3), 036401. (Journal details from the previous repository; the quote in D20 was read from the manuscript at personalpages.surrey.ac.uk/a.skeldon/DALEC_paper.pdf, dated 8 December 2014. Check against the published article.)
- Coleman, K., Prout, J. M. and Milne, A. E. (2024). RothC: a model for the turnover of carbon in soil. Model description. Rothamsted Research.
- Delahaies, S., Roulstone, I. and Nichols, N. (2017). Constraining DALEC v2 using multiple data streams and ecological constraints: analysis and application. Geoscientific Model Development 10, 2635-2650.
- DeLucia, E. H., Drake, J. E., Thomas, R. B. and Gonzalez-Meler, M. (2007). Forest carbon use efficiency: is respiration a constant fraction of gross primary production? Global Change Biology 13, 1157-1167.
- Douglas, et al. (2026). Regularisation of the 4DEnVar data assimilation method for calibration of land surface models. EGUsphere preprint.
- Famiglietti, C. A. et al. (2021). Optimal model complexity for terrestrial carbon cycle prediction. Biogeosciences 18, 2727.
- Fang, et al. Differentiable land model reveals global environmental controls on ecological parameters. arXiv:2411.09654.
- Fer, I. et al. (2018). Linking big models to big data: efficient ecosystem model calibration through Bayesian model emulation. Biogeosciences 15, 5801.
- Fox, A. et al. (2009). The REFLEX project: comparing different algorithms and implementations for the inversion of a terrestrial ecosystem model against eddy covariance data. Agricultural and Forest Meteorology 149, 1597-1615.
- Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808. (Sections used: 2.1 to 2.5, 4.1, 6.3, 7.1, 7.3, 10.5.)
- Hollinger, D. Y. and Richardson, A. D. (2005). Uncertainty in eddy covariance measurements and its application to physiological models. Tree Physiology 25.
- Ilvesniemi, H. et al. (2009). Long-term measurements of the carbon balance of a boreal Scots pine dominated forest ecosystem. Boreal Environment Research 14(4), 731-753. (Site carbon stocks and flows, Fig. 6 and Table 4; used for most priors.)
- Kallioinen, N., Paananen, T., Bürkner, P.-C. and Vehtari, A. (2024). Detecting and diagnosing prior and likelihood sensitivity with power-scaling. Statistics and Computing. doi:10.1007/s11222-023-10366-5 (volume and article number: check)
- Kattge, J. et al. (2011). TRY: a global database of plant traits. Global Change Biology 17, 2905-2935. (Cited in the official code for the lma range.)
- Kolari, P. (2010). Dissertationes Forestales 99, doi:10.14214/df.99.
- Kuppel, S. et al. (2012). Constraining a global ecosystem model with multi-site eddy-covariance data. Biogeosciences 9, 3757-3776.
- López-Blanco, E., Exbrayat, J.-F., Lund, M., Christensen, T. R., Tamstorf, M. P., Slevin, D., Hugelius, G., Bloom, A. A. and Williams, M. (2019). Evaluation of terrestrial pan-Arctic carbon cycling using a data-assimilation system. Earth System Dynamics 10, 233-255. (Uses the DALEC2 code released as Williams 2019.)
- Mahecha, M. D. et al. (2010). Global convergence in the temperature sensitivity of respiration at ecosystem level. Science 329, 838-840.
- Makela, A. et al. (2004). Acclimation of photosynthetic capacity in Scots pine to the annual cycle of temperature. Tree Physiology 24, 369-376. (Details to check.)
- Pastorello, G. et al. (2020). The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Scientific Data 7, 225.
- Pensa, M. and Jalkanen, R. (1999). Needle chronologies on Pinus sylvestris in northern Estonia and southern Finland. Silva Fennica 33(3), 171-177. doi:10.14214/sf.654 (pages to check)
- Pinnington, E. M. et al. (2016). Investigating the role of prior and observation error correlations in improving a model forecast of forest carbon balance using Four-dimensional Variational data assimilation. Agricultural and Forest Meteorology 228-229, 299-314.
- Pu, et al. (2025). MCI GPP. Scientific Data.
- Peltoniemi, M., Pulkkinen, M., Aurela, M., Pumpanen, J., Kolari, P. and Makela, A. (2015). A semi-empirical model of boreal-forest gross primary production, evapotranspiration, and soil water: calibration and sensitivity analysis. Boreal Environment Research 20, 151-171. (check)
- Richardson, A. D., Hollinger, D. Y., Burba, G. G. et al. (2006). A multi-site analysis of random error in tower-based measurements of carbon and energy fluxes. Agricultural and Forest Meteorology 136, 1-18. (check)
- Richardson, A. D. et al. (2010). Estimating parameters of a forest ecosystem C model with measurements of stocks and fluxes as joint constraints. Oecologia 164, 25-40.
- Schiestl-Aalto, P. (2017). Modelling intra- and inter-annual growth dynamics of Scots pine in the whole-tree carbon balance framework. Dissertationes Forestales 234.
- Schiestl-Aalto, P. et al. (2013). Duration of shoot elongation in Scots pine varies within the crown and between years. Annals of Botany. doi:10.1093/aob/mct180 (title and year to check)
- Smallman, T. L. and Williams, M. (2019). Description and validation of an intermediate complexity model for ecosystem photosynthesis and evapotranspiration: ACM-GPP-ETv1. Geoscientific Model Development 12, 2227-2253.
- Stettz, S. G. et al. (2022). Resolving temperature limitation on spring productivity in an evergreen conifer forest using a model-data fusion framework. Biogeosciences 19, 541-558.
- Suni, T. et al. (2003). Air temperature triggers the recovery of evergreen boreal forest photosynthesis in spring. Global Change Biology 9, 1410-1426. (Details to check.)
- Tramontana, G. et al. (2016). Predicting carbon dioxide and energy fluxes across global FLUXNET sites with regression algorithms. Biogeosciences 13, 4291.
- Tuomi, M. et al. (2009). Leaf litter decomposition: estimates of global variability based on Yasso07 model. Ecological Modelling 220, 3362-3371. (Cited via the previous repository; check.)
- Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing 27(5), 1413-1432.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. Bayesian Analysis 16(2), 667-718.
- Waring, R. H., Landsberg, J. J. and Williams, M. (1998). Net primary production of forests: a constant fraction of gross primary production? Tree Physiology 18, 129-134.
- Williams, M. et al. (1997). Predicting gross primary productivity in terrestrial ecosystems. Ecological Applications 7(3), 882-894. (Cited for the frost rule, D24; the rule itself not yet found in it.)
- Williams, M., Schwarz, P. A., Law, B. E., Irvine, J. and Kurpius, M. R. (2005). An improved analysis of forest carbon dynamics using data assimilation. Global Change Biology 11(1), 89-105. doi:10.1111/j.1365-2486.2004.00891.x (the evergreen DALEC; also the model reference of research_plan_v2.md).
- Williams, M. et al. (2009). Improving land surface models with FLUXNET data. Biogeosciences 6, 1341-1359.
- Williams, M. (2019). DALEC2 [software]. University of Edinburgh. https://doi.org/10.7488/ds/2504. (File `DALEC2 Feb19.f90`; read 2026-10-25, F33.)
- Yuan, Q. et al. (2025). Global carbon flux dataset generated by fusing remote sensing and multiple flux networks observation. Scientific Data.

Verify every citation against the original paper before using it in the thesis.
