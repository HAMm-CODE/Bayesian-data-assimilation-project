# Decisions

Every methodological and technical choice in this project, with its reason.
Each entry is marked with a status:

- **Decided**: settled and implemented.
- **Proposed**: preferred option chosen, not yet final.
- **Open**: still to be decided; the options are recorded.

Last updated: 2026-10-13.

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

**Size of the difference (real data, calibration block, 2026-10-12).** `t_mid - TA_F` averages -0.10 C; 90% of days lie between -1.05 and +0.78 C. Thesis sentence: the midpoint differs from the 24-hour mean by 0.1 C on average.

**Confirmed by the official code (2026-10-13).** The DALEC2 source code (CARDAMOM, `C/projects/CARDAMOM_MODELS/DALEC/DALEC_400/DALEC_400.c`, "Code used by Bloom et al., 2016. See also Bloom & Williams 2015") computes the decomposition multiplier as `exp(pars[9]*0.5*(maxtemp + mintemp))`, i.e. the midpoint. D20 matches the model as its authors run it, not only one paper's wording.

**Daily `TA_F` is still read**, but only for the D18 safety check (`drivers.temperature_extremes.daily_check_column`).

**Possible sensitivity check.** Because the two readings are both defensible, a short run with `TA_F` instead of `t_mid` would show whether the choice moves the posterior. Worth doing only if Figure H shows a sizeable difference. **Result: not needed.** An average difference of 0.1 C is small next to the 6.3 C daily range and the seasonal swing, so the run is dropped.

---

## Part 2b. Model structure (decided 2026-10-13; details in MODEL.md)

**Reference used.** Besides the papers, the official DALEC2 code (CARDAMOM, `DALEC_400`) was read: `DALEC_400.c`, `PARS_INFO_400.c`, `DALEC_ALL/ACM.c`, `DALEC_ALL/offset.c` (github.com/CARDAMOM-framework/CARDAMOM). Its header says it is the code used by Bloom et al. (2016), with Bloom and Williams (2015) as the model description. Where paper and code disagree, the code shows what the authors actually run. Caveat: its prior ranges are the 2016 ones and differ from the 2015 Table 1 (relevant to priors, not structure).

### D21. ACM in the Chuter et al. (2015) form (M1)
Day length from latitude and solar declination, conductance `|psi_d|^a10 / (0.5 Tr + a6 R_tot)`, temperature on `Tmax`. **Why:** `ACM.c` in the official code computes exactly this. The Williams et al. (1997) form has no latitude in its day length and would overstate winter GPP at 62 N. This was an assumption in the previous repository; it is now confirmed.

### D22. ACM coefficients are DALEC2's fixed constants, at full precision (M2)
`a2 = 0.0156935, a3 = 4.22273, a4 = 208.868, a5 = 0.0453194, a6 = 0.37836, a7 = 7.19298, a8 = 0.011136, a9 = 2.1001, a10 = 0.789798, psi_d = -2.0 MPa, R_tot = 1.0` (`DALEC_400.c`, `constants[]` and `gpppars[]`). **Why:** these are not a site choice; DALEC2 fixes them. They equal the "Loobos" set of Chuter et al. (2015), which the previous repository used in rounded form (for example 0.0156). They were fitted outside the boreal zone; stated limitation.

### D23. Solar day counted from 21 December (`doy + 10`) (M3)
**Why:** the declination in `ACM.c` uses `(yearday + 10)`, matching Chuter et al. Section II.A. Without it the seasonal cycle of day length is 10 days late.

### D24. No frost cut-off in the main model; kept as a switch (M4)
The previous repository stopped GPP when `(Tmax + Tmin)/2 < -2 C`. **Why not by default:** no such rule exists in the official code, and the cited precedent (Williams et al. 1997, Oregon conifers) has not been verified. Without it the thesis can state that DALEC2 is calibrated as distributed. (Note: the code's `gpppars[8] = -2.0` is the leaf water potential `psi_d`, not a frost threshold.) **Kept as a switch** because winter days carry most of the likelihood weight and DALEC2 lacks cold limitation on photosynthesis (F9, Stettz et al. 2022); it is tested once in the prior predictive check and once as a sensitivity run.

### D25. A8 corrected: leaf fall anchored on `d_fall` (M5)
`Phi_fall` uses `t - d_fall + psi_f`, not the printed `t - cr_fall + psi_f`. **Why:** the official code anchors leaf fall on its "Fday" parameter (`pars[14]`), the leaf-fall day, mirroring how leaf onset is anchored on "Bday". The printed A8 is a typo. This also removes the previous repository's inert `d_fall`.

### D26. All five phenology parameters are fitted, with site-based ranges (M6)
`d_onset`, `cr_onset`, `d_fall`, `cr_fall`, `c_lf` are sampled. **Why:** RQ1 and RQ3 ask what NEE can constrain; fixing parameters would hide the answer. The ranges are set in step 3 from Scots pine knowledge at Hyytiala, not the 1 to 365 day published range.

### D27. Starting pools are sampled, with priors centred on site measurements (M7)
All six `C(0)` are parameters with priors from Hyytiala stock measurements (Ilvesniemi et al. 2009; set in step 3). **Why:** the simplest first model (Gelman et al. 2020, Section 2.1). **Fallback:** the previous repository's reparameterisation (sample soil respiration and its split, derive the litter and soil pools) if the sampler struggles with the strong link between a turnover rate and its pool.

### D28. Allocation through a Dirichlet split (M8)
`f_auto` plus a Dirichlet split of the remaining `1 - f_auto` across labile, foliage, roots and wood. **Why:** the official code instead splits step by step (foliage share of NPP, then labile share of the rest, then roots, wood taking the remainder). Both forms reach exactly the same set of allocations, so the model is the same; only the way the prior is written differs. The Dirichlet makes site-based priors simple and never produces negative fractions. This also re-frames correction 1 in `MODEL.md`: the paper's footnote is a shorthand of the code's step-by-step split.

### D29. LAI convention moved to the priors (M9)
It only changes the prior range of `lma`, so it is decided with that prior in step 3.

### D30. Phenology uses a continuous day counter (M10)
The leaf onset and leaf fall pulses use `t`, the number of days since the start of the run, not the day of year. **Why:** the official code does this (`MET[:,0]`, "projday"). With a 365.25-day cycle, a day-of-year counter that restarts every 1 January causes a small jump each year. Day length (B5) still uses the day of year, as in the code.

### D31. Leaf-fall offset `psi` from the official code's polynomial (A9)
`psi` is computed with the sixth-order polynomial in `log(L - 1)`, `L = 1 / c_lf`, from `DALEC_ALL/offset.c`. **How it was settled:** the polynomial matches the equation `2 sqrt(pi) psi = log(1 - c_lf) exp(-psi^2)` within 0.004 for leaf lifespans 1.001 to 8 years; it does not match the previous repository's form (`log(1 - c_lf)` on the other side), which gave an offset about 8 times too large at `c_lf = 0.25` (about 24 days too early or late for a 60-day fall period). The A7 constants 0.6245 and 6.9088 come out of the same polynomial at `L = 1.001`, which confirms the reading. **Why the polynomial, not an exact root:** `c_lf` is now sampled (D26), so `psi` must change smoothly with it for NUTS gradients; the polynomial is differentiable, the root solve is not without extra work, and it is what DALEC2's authors run. Limit: valid for `L` in 1.001 to 8 years only, which covers Scots pine needles (3 to 5 years). **Still to check:** the printed A9 in the PDF, to tell a paper typo from a copying error.

### D32. Parameter scales (Gelman et al. 2020, Section 2.3)
Each parameter's sampling scale follows from its range by one rule: **log** if upper / lower >= 10, **interval** (logit) if less or if the parameter is a day-of-year position, **simplex** for the allocation split. Result: all five turnover rates, `cr_onset`, `lma`, `ceff` and all six starting pools on the log scale; `f_auto`, the temperature exponent, `d_onset`, `d_fall`, `cr_fall` and `c_lf` on the interval scale; the four allocation shares on the simplex. **Why:** NUTS moves best when every parameter is unbounded and changes on a similar scale; the published ranges span from a factor of 4 to a factor of 10,000. The papers give ranges, not scales, so the scale is a modelling choice; using one rule keeps it simple and explainable. Re-applied to the final ranges after step 3. In PyMC most of this is automatic: a log-normal prior is sampled on the log scale, a bounded prior on the interval scale.

### D33. Prior strategy: smooth, weakly informative priors (Gelman et al. 2020, Section 2; step 3)
Every prior is a smooth distribution on the parameter's scale (D32), centred on a site measurement or literature value and wide enough that the published range stays possible but its extremes are unlikely. No flat priors with hard edges. **Why:** hard edges are walls for NUTS; the previous repository's unusable run had eight parameters pinned at their prior bounds. Smooth priors also state the belief honestly: values near the measurement are likely, extremes possible but rare. **Options not taken:** flat priors on narrowed ranges (the previous repository's approach), and a mix of flat and smooth. **Cost:** each prior needs a stated centre and width with a source; recorded parameter by parameter below.

### D34. Priors, group 1: allocation
- `f_auto ~ Beta(25.5, 24.5)`: centre 0.51, sd 0.07, 95% 0.37 to 0.65; under 0.4% outside the published 0.3 to 0.7. **Source:** Hyytiala tree growth 505 of GPP about 1,030 g C m-2 yr-1 (Ilvesniemi et al. 2009, Fig. 6) gives 0.51, site range 0.37 to 0.60; Waring et al. (1998) NPP/GPP 0.47 +/- 0.04 (f_auto 0.53); DeLucia et al. (2007) NPP/GPP 0.23 to 0.83 across forests, lowest in boreal forests (f_auto above 0.5).
- `(f_lab, f_fol, f_roo, f_woo) = (1 - f_auto) x Dirichlet(13.06, 13.06, 15.26, 44.27)`: shares of growth 0.152, 0.152, 0.178, 0.517 from the measured needle litterfall 154 (split equally between labile and foliage), below-ground litter 90 and wood growth 261 g C m-2 yr-1 (Ilvesniemi et al. 2009, Fig. 6); total concentration 85.65 derived from the measured spread of litterfall and GPP (previous repository, section 11). As shares of GPP: medians 0.073, 0.073, 0.085, 0.251.
- **Why this informative:** the previous repository's flat split sent too much carbon to the needles and LAI ran away to 48 (site value about 3); Gelman et al. (2020, Section 7.3): "as the model becomes more complex, all of the priors need to become tighter".
- **Caveats:** the equal labile/foliage split is a judgement (the data constrain only their sum); tower GPP includes the understory but the growth flows are trees only.
- **Planned check (Gelman et al. 2020, Section 6.3):** after fitting, compare prior and posterior spread; if the posterior of `f_auto` is about as wide as its prior, refit once with sd 0.10 as a sensitivity run.

### D35. Priors, group 2: turnover rates
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
1. **Judge the priors together, not one by one** (Sections 2.4, 7.3: independent priors "may jointly be strongly informative"). The prior predictive check (`scripts/03_prior_predictive.py`) will run DALEC2 with draws from all priors at once. Its pass criteria are written down **before** it is run (annual GPP, ecosystem and soil respiration, NEE, peak LAI, and pool sizes after 14 years), so the check cannot be tuned after seeing it.
2. **Watch the soil carbon used twice:** in `theta_som` here and in the starting soil pool (D27). The check above shows whether the implied soil respiration is sensible; if not, the previous repository's reparameterisation (D27 fallback) is used.
3. **After fitting, measure how much each prior mattered** (Section 6.3): compare prior and posterior spread for every parameter. Parameters whose posterior is about as wide as their prior were not informed by NEE; this is a direct result for RQ1 (expected for the slow `theta_woo` and `theta_som`). Only those get a sensitivity refit with a wider prior.

### D36. Priors, group 3: temperature exponent (confirmed 2026-10-13)
`Theta ~ LogNormal(median 0.0336, sdlog 0.35)`: 95% 0.017 to 0.067, i.e. Q10 = exp(10 Theta) of 1.18 to 1.95 (median Q10 1.40). 3.7% below and 0.7% above the published 0.018 to 0.08. **Source:** Mahecha et al. (2010): ecosystem respiration across 60 FLUXNET sites converges to Q10 = 1.4 +/- 0.1 against air temperature once seasonal confounding is removed; the previous repository's estimate from Hyytiala winter fluxes was Theta = 0.0366 (Q10 1.44), close to the median. **Why wider than Mahecha's +/- 0.1:** DALEC2's single exponent also absorbs effects the model lacks (soil frost, moisture), so a site-level value may differ from the global one. **Scale changed from interval to log** (D32 allows this after step 3): the prior is now positive-only and smooth without needing the published bounds, which it already respects on 96% of draws.

### D37. Priors, group 4: leaf timing (phenology)
In DALEC2 terms at an evergreen pine site, the onset pulse is the **spring flush of new shoots and needles** and the fall pulse is the **autumn drop of the oldest needles**. Note: this is **not** the NEE spring onset of Figure G (median day 100), which is the recovery of photosynthesis in needles that already exist; DALEC2 cannot represent that recovery (F9), and its leaf onset should not be tuned to it.

| Parameter | Prior | 95% | Basis |
|---|---|---|---|
| `d_onset` | Normal(145, 15) | day 116 to 174 | Scots pine shoot elongation in southern Finland starts in about mid to late May; with `cr_onset` 40 the release peaks around day 170 (mid-June), during needle elongation. **To confirm** with Schiestl-Aalto et al. (2013, Annals of Botany) or the dissertation (Schiestl-Aalto 2017), which give timing in heat sums |
| `cr_onset` | LogNormal(median 40 d, sdlog 0.35) | 20 to 79 days | length of shoot and needle elongation, roughly May to July; inside the published 10 to 100 |
| `d_fall` | Normal(265, 15) | day 236 to 294 | day of maximum needle fall (D25: the paper says the offset places the maximum loss rate on `d_fall`); Scots pine needle fall in Finland mostly late August to October. **To confirm** from monthly litterfall data |
| `cr_fall` | LogNormal(median 60 d, sdlog 0.35) | 30 to 119 days | spread of needle fall over about two months; inside the published 20 to 150 |
| `c_lf` | Beta(14.51, 39.24) | 0.16 to 0.40 (needle life 2.5 to 6.2 years) | mean needle retention 3.4 to 4.2 years in southern Finland (Pensa and Jalkanen 1999); mean 0.27 (3.7 years). 0.3% below 0.125, the edge of the D31 polynomial |

Day positions use a plain normal prior (smooth, no bounds needed because the window is far from day 0 and day 365); `cr_fall` moves from interval to log scale (D32 re-applied). Days are counted on the continuous counter (D30), with day 1 = 1 January 1997.

### D38. Priors, group 5: canopy (`lma`, `ceff`), and the LAI convention (closes D29)
- **`lma ~ LogNormal(median 185 g C m-2, sdlog 0.25)`**, 95% 113 to 302. **Basis:** foliage carbon from needle litterfall / `c_lf` = 154 / 0.26, about 590 g C m-2 (Ilvesniemi et al. 2009), divided by site LAI. Kolari (2010) gives all-sided LAI 6.5 to 8; as hemisurface LAI (divide by 2) this is 3.2 to 4.0, giving `lma` about 165; as projected LAI (divide by 2.57) 2.5 to 3.1, giving about 210. **The LAI convention question (M9) is resolved by putting it inside the prior:** the median sits between the two conventions and the 95% range covers both. One run instead of two sensitivity runs, and the uncertainty is stated honestly (Gelman et al. 2020, Section 7.3: uncertainty belongs in the prior). Note: about a third of the prior lies above the official code's upper limit of 200 (from Kattge et al. 2011), but inside the 2015 paper's 400.
- **`ceff ~ LogNormal(median 16.5, sdlog 0.41)`**, 95% 7.4 to 37. **Basis:** `ceff` replaces `p11 x N` (nitrogen-use efficiency x foliar nitrogen per leaf area). `p11 = 7.4` belongs to the same ACM fit as the fixed coefficients (D22, Chuter et al. 2015, Loobos); with needle N between about 1 and 5 g N m-2 leaf area this gives 7.4 to 37, which also contains the Loobos value 29.6 and most of the previous repository's REFLEX range 5 to 20. **To check:** a Hyytiala needle N per area value would narrow this.

### D39. Priors, group 6: starting pools on 1 January 1997 (implements D27)
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

### Summary of step 3
24 sampled quantities: `f_auto`, the 4-share Dirichlet (3 free), 4 turnover rates, `h`, `Theta`, 5 phenology, `lma`, `ceff`, 6 starting pools. Every prior is smooth (D33) and has a stated source or reason; three values are marked **to confirm** (phenology timing, needle fall timing, needle N). Next: judge them together in the prior predictive check (D35, way forward 1).

### D40. Prior predictive check: pass criteria, written before running (Gelman et al. 2020, Section 2.4)
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

**F10. First real-data run of the driver figures (2026-10-12).**
- The half-hourly safety check (D18) passed: every day has 48 valid half-hours, and the half-hourly mean matches the daily `TA_F` within 0.05 C on every day.
- True daily temperature range averages 6.3 C over 1997-2010. All model inputs are complete in the calibration block (no missing days).
- Spring onset with the 10-day rule (D16): median day 100, range 66 to 137 over the calibration years; evaluation years 93 to 111. See OPEN-3.

**F11. Days are not independent (Figure K).** NEE anomalies are correlated with the next day (+0.37), two days later (+0.22) and a week later (+0.08); they reach zero by about two weeks. If days were independent, values would stay within +/-0.03. The Gaussian likelihood (L2) assumes independent days, so it treats the data as more informative than it is and the posterior may be too narrow. Recorded as a limitation; an error model with day-to-day correlation (for example AR(1), where each day's error carries part of the previous day's) is a possible later sensitivity analysis.

**F6. Evaluation block.** 608 of 1,461 evaluation days (41.6%) are excluded only because RANDUNC is missing: all of 2014 and the first half of 2011.

---

## Part 5. Open items

**OPEN-1. Evaluation mask.** Evaluation needs sigma only for probability-based checks (for example, predictive interval coverage). Error measures like RMSE need only NEE and QC. Using QC alone keeps 1,418 instead of 810 evaluation days (adds 2014 and early 2011). Decide when designing the evaluation chapter.

**OPEN-3. Spring onset number.** The previous EDA found a median observed onset of day 96 (range 58 to 117) with the 10-day rule, while project notes mention an observed onset of day 110 in the model comparison. The new code gives a median of day 100 (range 66 to 137) for 1997-2010 (F10), so three numbers now exist. These probably come from different definitions or data. Pick one definition (D16 uses the 10-day rule) and use it consistently in the thesis.

**OPEN-4. Daily temperature extremes (F7).** Closed 2026-10-12: derived from the half-hourly file (D18).

**OPEN-5. Model structure.** Closed 2026-10-13: M1 to M8 and M10 decided (D21 to D28, D30); M9 moved to the priors (D29).

### Carried over from the previous repository, not yet re-assessed

These will be reviewed when the work reaches them (allocation, ACM form, LAI convention and `d_fall` were re-assessed on 2026-10-13: D21, D25, D28, D29):

- No EDC constraints (they break differentiability for NUTS).
- No VPD driver.
- NUTS as the only sampler.

---

## References

- Bacour, C. et al. (2023). Assimilation of multiple datasets results in large differences in regional- to global-scale NEE and GPP budgets simulated by a terrestrial biosphere model. Biogeosciences 20, 1089-1111.
- Bloom, A. A. and Williams, M. (2015). Constraining ecosystem carbon dynamics in a data-limited world: integrating ecological "common sense" in a model-data fusion framework. Biogeosciences 12, 1299-1315.
- Bloom, A. A., Exbrayat, J.-F., van der Velde, I. R., Feng, L. and Williams, M. (2016). The decadal state of the terrestrial carbon cycle: global retrievals of terrestrial carbon allocation, pools, and residence times. PNAS 113(5), 1285-1290. (Version of DALEC2 in the CARDAMOM code `DALEC_400`.)
- CARDAMOM framework, source code. github.com/CARDAMOM-framework/CARDAMOM, files `C/projects/CARDAMOM_MODELS/DALEC/DALEC_400/DALEC_400.c`, `PARS_INFO_400.c`, `DALEC_ALL/ACM.c`, `DALEC_ALL/offset.c` (read 2026-10-13).
- Brunner, I. et al. (2013). Fine-root turnover rates of European forests revisited: an analysis of data from sequential coring and ingrowth cores. Plant and Soil 362, 357-372.
- Coleman, K., Prout, J. M. and Milne, A. E. (2024). RothC: a model for the turnover of carbon in soil. Model description. Rothamsted Research.
- Chuter, A. M., Aston, P. J., Skeldon, A. C. and Roulstone, I. (2015). A dynamical systems analysis of the data assimilation linked ecosystem carbon (DALEC) models. Chaos 25(3), 036401. (Journal details from the previous repository; the quote in D20 was read from the manuscript at personalpages.surrey.ac.uk/a.skeldon/DALEC_paper.pdf, dated 8 December 2014. Check against the published article.)
- Delahaies, S., Roulstone, I. and Nichols, N. (2017). Constraining DALEC v2 using multiple data streams and ecological constraints: analysis and application. Geoscientific Model Development 10, 2635-2650.
- DeLucia, E. H., Drake, J. E., Thomas, R. B. and Gonzalez-Meler, M. (2007). Forest carbon use efficiency: is respiration a constant fraction of gross primary production? Global Change Biology 13, 1157-1167.
- Douglas, et al. (2026). Regularisation of the 4DEnVar data assimilation method for calibration of land surface models. EGUsphere preprint.
- Famiglietti, C. A. et al. (2021). Optimal model complexity for terrestrial carbon cycle prediction. Biogeosciences 18, 2727.
- Fang, et al. Differentiable land model reveals global environmental controls on ecological parameters. arXiv:2411.09654.
- Fer, I. et al. (2018). Linking big models to big data: efficient ecosystem model calibration through Bayesian model emulation. Biogeosciences 15, 5801.
- Fox, A. et al. (2009). The REFLEX project: comparing different algorithms and implementations for the inversion of a terrestrial ecosystem model against eddy covariance data. Agricultural and Forest Meteorology 149, 1597-1615.
- Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808 (Section 10.5, pp. 53-54; Section 7.1).
- Hollinger, D. Y. and Richardson, A. D. (2005). Uncertainty in eddy covariance measurements and its application to physiological models. Tree Physiology 25.
- Kattge, J. et al. (2011). TRY: a global database of plant traits. Global Change Biology 17, 2905-2935. (Cited in the official code for the lma range.)
- Kolari, P. (2010). Dissertationes Forestales 99, doi:10.14214/df.99.
- Kuppel, S. et al. (2012). Constraining a global ecosystem model with multi-site eddy-covariance data. Biogeosciences 9, 3757-3776.
- Mahecha, M. D. et al. (2010). Global convergence in the temperature sensitivity of respiration at ecosystem level. Science 329, 838-840.
- Pensa, M. and Jalkanen, R. (1999). Needle chronologies on Pinus sylvestris in northern Estonia and southern Finland. Silva Fennica 33(3), 171-177. doi:10.14214/sf.654 (pages to check)
- Pastorello, G. et al. (2020). The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Scientific Data 7, 225.
- Pu, et al. (2025). MCI GPP. Scientific Data.
- Pinnington, E. M. et al. (2016). Investigating the role of prior and observation error correlations in improving a model forecast of forest carbon balance using Four-dimensional Variational data assimilation. Agricultural and Forest Meteorology 228-229, 299-314.
- Richardson, A. D. et al. (2010). Estimating parameters of a forest ecosystem C model with measurements of stocks and fluxes as joint constraints. Oecologia 164, 25-40.
- Smallman, T. L. and Williams, M. (2019). Description and validation of an intermediate complexity model for ecosystem photosynthesis and evapotranspiration: ACM-GPP-ETv1. Geoscientific Model Development 12, 2227-2253.
- Schiestl-Aalto, P. et al. (2013). Duration of shoot elongation in Scots pine varies within the crown and between years. Annals of Botany. doi:10.1093/aob/mct180 (title and year to check)
- Schiestl-Aalto, P. (2017). Modelling intra- and inter-annual growth dynamics of Scots pine in the whole-tree carbon balance framework. Dissertationes Forestales 234.
- Stettz, S. G. et al. (2022). Resolving temperature limitation on spring productivity in an evergreen conifer forest using a model-data fusion framework. Biogeosciences 19, 541-558.
- Tuomi, M. et al. (2009). Leaf litter decomposition: estimates of global variability based on Yasso07 model. Ecological Modelling 220, 3362-3371. (Cited via the previous repository; check.)
- Waring, R. H., Landsberg, J. J. and Williams, M. (1998). Net primary production of forests: a constant fraction of gross primary production? Tree Physiology 18, 129-134.
- Tramontana, G. et al. (2016). Predicting carbon dioxide and energy fluxes across global FLUXNET sites with regression algorithms. Biogeosciences 13, 4291.
- Williams, M. et al. (2009). Improving land surface models with FLUXNET data. Biogeosciences 6, 1341-1359.
- Yuan, Q. et al. (2025). Global carbon flux dataset generated by fusing remote sensing and multiple flux networks observation. Scientific Data.

Verify every citation against the original paper before using it in the thesis.