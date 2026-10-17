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

### D57. Priors of the cold-limit model: `ceff` returns to its first prior, D38 (Decided as a test, 2026-10-16)
**Why.** F22: with `ceff` near 10, ACM photosynthesis hardly rises with sunlight, and the fit escapes to an implausible canopy. D42 had lowered `ceff` from D38's median 16.5 (95% 7.4 to 37, built from p11 x needle N before any data were seen) to median 10 (95% 5 to 20). Its only stated reason: with no cold limit, higher values make annual GPP implausibly high, because the canopy photosynthesises through the winter. The cold-limit model removes that reason, so the narrowing has to be checked again for that model.

**What.** In the cold-limit model only, `ceff ~ LogNormal(median 16.5, sdlog 0.41)` (D38). Everything else as before. The main model keeps D42.
- Code: `COLD_MODEL_PRIOR_CHANGES` and `priors_for(cold)` in `priors.py`; `build_model(..., cold_limit=True)` uses them. Each fit now saves its prior list with the results (`priors` note in the .nc file), so the parameter table and `diagnose_chains.py` always use the priors the fit was made with. Fits saved earlier carry no note; they all used the main priors.

**Check before any fit (prior predictive, script 03 `--cold`, no tower NEE used).** Same 1,000-draw design and D40 checks as F13 to F15, with the soft foliage constraint (D44) as fitted. **Rule, fixed before running:** use the D38 `ceff` prior in the cold model if Tier 2 is at least 50% and Tier 1 at most 7% (the main model was accepted at 6.4%, F15; with about 800 effective draws the Monte Carlo error of such a rate is about 0.9 percentage points). Otherwise keep D42 for the cold model too (option B of F22). The script also runs the same draws with the D42 `ceff` prior, for comparison.

**Then:** a quick fit of the cold model (`--cold`), judged like F21 and F22. If it is believable, the cold model goes through the fake-data check and the full fit (D55).

**Gelman et al. (2020) on this step:**
- Section 2.4: prior predictive checks "provide a way to refine the model without using the data multiple times". That is why the decision rests on the prior predictive check and on the D38 and D42 sources, not on the fit.
- Section 7.3: "Priors must be specified for each model within a workflow." The cold-limit model is a new model, so its priors are set and checked on their own.
- Section 7.3 also notes that "often as the model becomes more complex, all of the priors need to become tighter." Here one prior is loosened. The reason is specific: the tight `ceff` prior did the job of the missing cold limit, and the cold limit now does that job itself.
- Section 6.3: sensitivity analysis by "refitting the model with multiple priors". The cold fit with the D42 prior (F21) is kept and reported as that comparison.
- **Stated openly:** this change was considered after seeing that the fit prefers high `ceff` (F21, F22). It is not a free choice of the best-fitting prior: it returns to a prior written down before any data were seen, for a reason recorded in D42 at the time, and it must pass a check that uses no tower data.

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
| `f_fol` | not printed then (0.015 to 0.017 in S1 and S4, F20) | 0.071 (prior median about 0.07) |
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

## Part 5. Open items

### Still open

**OPEN-1. Evaluation mask.** Evaluation needs sigma only for probability-based checks (for example, predictive interval coverage). Error measures like RMSE need only NEE and QC. Using QC alone keeps 1,418 instead of 810 evaluation days (adds 2014 and early 2011). Decide when designing the evaluation chapter.

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

**OPEN-8. How to respond to F17 (real NEE pulls the parameters out of the believable region).** Options, all already recorded as planned changes, none tuned to the result: (a) add model error to sigma (L1 option C, Gelman et al. 2020, Sections 7.1 and 10.5), so the likelihood stops treating DALEC2 as able to match each day within RANDUNC; (b) the frost switch (D24), aimed at the spring misfit; (c) both; (d) a Laplace error (L2). Plan: quick fits of (a), (b) and (c) first, each judged with `diagnose_chains.py` (are the chains believable and do they agree?). **Status 2026-10-16 (F18):** (a) removes the separate modes and is now the main error model (L1, changed 2026-10-16); (b) alone does nothing. The remaining implausible canopy leaves three routes (Gelman et al. 2020, Sections 5.2, 6.1, 7): tolerate and report it; simplify (fix the canopy parameters at site values); or expand (a cold limit on photosynthesis, Stettz et al. 2022). To discuss with the supervisor. Gelman et al. (2020, Section 5.9) warn against "haphazardly tweaking the model until computational problems disappear", so every option must have its own reason, and all results are reported, not only the one that works.

### Planned sensitivity runs (decided in principle, to run after the main fit)

| Run | What changes | Entry |
|---|---|---|
| Frost switch | GPP stopped when the midpoint temperature is below -2 C | D24 |
| RANDUNC only (the original main model) | sigma_t = RANDUNC_t, no model error | L1 (main model changed 2026-10-16) |
| Wider `f_auto` prior | sd 0.10 instead of 0.07, only if the posterior is about as wide as the prior | D34 |
| Wider priors where NEE gave no information | only for parameters whose posterior is about as wide as the prior | D35, way forward 3 |
| Day-to-day correlated errors | for example AR(1) errors, if time allows | F11 |

### Closed

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
- Kattge, J. et al. (2011). TRY: a global database of plant traits. Global Change Biology 17, 2905-2935. (Cited in the official code for the lma range.)
- Kolari, P. (2010). Dissertationes Forestales 99, doi:10.14214/df.99.
- Kuppel, S. et al. (2012). Constraining a global ecosystem model with multi-site eddy-covariance data. Biogeosciences 9, 3757-3776.
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
- Yuan, Q. et al. (2025). Global carbon flux dataset generated by fusing remote sensing and multiple flux networks observation. Scientific Data.

Verify every citation against the original paper before using it in the thesis.