# The DALEC2 model

What the model is, every equation it runs, and every parameter it has. It is the
starting point for the priors and the code.

This is **step 1** of Gelman et al. (2020), Section 2 ("Before fitting a model"):
write the model down before fitting it. Steps 2 and 3, choosing the scale and
prior for each parameter, fill the empty columns of the parameter table in
Section 6. Open model choices are listed in Section 8 and are moved to
`DECISIONS.md` once decided.

Status: draft, 2026-10-12. Equations transcribed from Bloom and Williams (2015,
Appendix A) and Chuter et al. (2015, Appendix B), cross-checked against the
previous repository. Every equation must still be checked against the PDFs
before it goes into the thesis (see Section 9).

---

## 1. What DALEC2 does, in plain words

DALEC2 is a "box model" of the carbon in a forest. It keeps track of six
**pools** (stores of carbon, in g C m-2), and moves carbon between them once per
day:

```
                      sunlight, temperature, CO2
                                 |
                                 v
                        GPP (photosynthesis)
                                 |
        +-------------+----------+-----------+------------+
        |             |          |           |            |
   f_auto (back   f_lab      f_fol       f_roo        f_woo
   to the air)      |          |           |            |
        |        LABILE ---> FOLIAGE      ROOTS        WOOD
        |        (stored    (needles)       |            |
        |         sugar)       |            |            |
        |                      +---> LITTER <+            |
        |                              |                  |
        |                              +---> SOIL <-------+
        |                              |      |
        v                              v      v
       Ra                       Rh (decomposition, faster when warm)
```

- **GPP** (gross primary production): carbon the trees take from the air by
  photosynthesis each day.
- **Allocation fractions** (`f_...`): how each day's GPP is shared out. One part
  (`f_auto`) is breathed straight back out by the trees; the rest grows the four
  plant pools.
- **Turnover rates** (`theta_...`): the fraction of a pool that leaves it each
  day (roots and wood die and become litter or soil; litter and soil decompose).
- **Respiration**: carbon going back to the air. `Ra` comes from the trees
  (autotrophic), `Rh` from decomposing litter and soil (heterotrophic). Their
  sum is ecosystem respiration, `Reco`.
- **NEE** (net ecosystem exchange) = `Reco - GPP`. Negative NEE means the forest
  took up carbon that day. NEE is what the tower measures and what the model is
  fitted to.

The model is **deterministic**: the same parameters and inputs always give the
same NEE. Randomness enters only through the likelihood (`DECISIONS.md` L1, L2).

## 2. Inputs and fixed constants

**Daily inputs** (`DECISIONS.md` D18 to D20):

| Symbol | Meaning | Unit | Column |
|---|---|---|---|
| `Tmax` | daily highest air temperature | degC | from half-hourly `TA_F` |
| `Tmin` | daily lowest air temperature | degC | from half-hourly `TA_F` |
| `T` | respiration temperature, `(Tmax + Tmin) / 2` | degC | computed |
| `Tr` | daily temperature range, `Tmax - Tmin` | degC | computed |
| `I` | incoming shortwave radiation, daily total | MJ m-2 d-1 | `SW_IN_F` x 0.0864 |
| `Ca` | atmospheric CO2 | umol mol-1 | `CO2_F_MDS` |
| `doy` | day of year | day | from the date |

**Fixed site constant:** latitude 61.85 N (sets day length).

**Fixed model constants** (not fitted): the ACM coefficients `a2` to `a10`,
`psi_d` and `R_tot` (Section 7), and the phenology constants `s`, 6.9088 and
0.6245 (Section 3.3).

## 3. Equations

Time steps are days, `t` to `t+1`. Every right-hand side uses the pools at day
`t`, and all six pools are updated **together** (updating one pool first and
then using its new value in another equation would give wrong numbers).

### 3.1 Carbon pools (Bloom and Williams 2015, A1 to A6)

```
(A1) labile   C_lab(t+1) = (1 - Phi_onset) C_lab + f_lab GPP
(A2) foliage  C_fol(t+1) = (1 - Phi_fall) C_fol + Phi_onset C_lab + f_fol GPP
(A3) roots    C_roo(t+1) = (1 - theta_roo) C_roo + f_roo GPP
(A4) wood     C_woo(t+1) = (1 - theta_woo) C_woo + f_woo GPP
(A5) litter   C_lit(t+1) = [1 - (theta_lit + theta_min) e^(Theta T)] C_lit
                           + theta_roo C_roo + Phi_fall C_fol
(A6) soil     C_som(t+1) = [1 - theta_som e^(Theta T)] C_som
                           + theta_woo C_woo + theta_min e^(Theta T) C_lit
```

`e^(Theta T)` is the **temperature multiplier**: decomposition speeds up
exponentially as it gets warmer. At `T = 0 degC` it equals 1, so the `theta`
rates of litter and soil are their rates at 0 degC. `T` is the midpoint
temperature (D20).

`theta_min` moves carbon from litter into soil. It is a transfer, not a loss, so
it does not appear in respiration.

### 3.2 Fluxes (A10 to A12)

```
(A10) Reco = f_auto GPP + (theta_lit C_lit + theta_som C_som) e^(Theta T)
             \_ Ra _/     \______________ Rh _______________/
(A11) NEE  = Reco - GPP
(A12) LAI  = C_fol / c_lma
```

`LAI` (leaf area index) is the leaf area per ground area. It is not an input: it
comes from the foliage pool, so the model grows and loses its own canopy.

**Built-in check (carbon balance).** Adding A1 to A6, every transfer between
pools cancels, which leaves

```
change in total carbon = GPP - Reco = -NEE      (exactly, every day)
```

The code will test this on every run.

### 3.3 Leaf phenology (A7 to A9)

Phenology is the timing of leaf growth and leaf fall. Each day, a fraction
`Phi_onset` of the labile pool moves into the foliage, and a fraction `Phi_fall`
of the foliage drops into the litter. Both are bell-shaped pulses once a year.

```
(A7) Phi_onset = sqrt(2/pi) (6.9088 / cr_onset)
                 exp( -[ sin((doy - d_onset - 0.6245 cr_onset) / s) sqrt(2) s / cr_onset ]^2 )

(A8) Phi_fall  = sqrt(2/pi) (-log(1 - c_lf) / cr_fall)
                 exp( -[ sin((doy - cr_fall + psi_f) / s) sqrt(2) s / cr_fall ]^2 )
                 (as printed; see correction 2 in Section 4)

     s = 365.25 / pi

(A9) psi solves  2 sqrt(pi) log(1 - c_lf) psi - exp(-psi^2) = 0,   psi_f = psi cr_fall / sqrt(2)
```

`psi` is not a parameter: it is calculated from `c_lf`. The previous repository
solved A9 numerically. **Check against the PDF:** one text extraction of the
paper showed `exp(-psi^2 / 2)` in A9 instead of `exp(-psi^2)`.

Hyytiala is an **evergreen** pine forest: it keeps its needles all year and
sheds only the oldest ones. How these deciduous-style pulses should behave there
is an open choice (M6).

### 3.4 Photosynthesis: the Aggregated Canopy Model (ACM)

Bloom and Williams (2015) give GPP only as a function from ACM (Williams et al.
1997). The form DALEC actually uses is written out by Chuter et al. (2015,
Appendix B, Eq. 12), which differs from the 1997 paper, most importantly in
using true day length from latitude. The previous repository used Chuter's form
(M1).

In solve order:

```
(B1) canopy conductance   g_c = |psi_d|^a10 / (0.5 Tr + a6 R_tot)
(B3) photosynthate        p   = ceff LAI e^(a8 Tmax) / g_c
(B2) internal CO2         q   = a3 - a4
                          C_i = 0.5 [ Ca + q - p + sqrt( (Ca + q - p)^2 - 4 (Ca q - a3 p) ) ]
     CO2-limited rate     p_d = g_c (Ca - C_i)
(B4) light-use efficiency E_0 = a7 LAI^2 / (LAI^2 + a9)
     light-limited rate   p_I = E_0 I p_d / (E_0 I + p_d)
(B5) solar declination    delta = -0.408 cos( 2 pi (doy + 10) / 365 )
(B6) day length (hours)   h = 24 arccos( -tan(latitude) tan(delta) ) / pi
(Eq. 12) GPP = p_I (a2 h + a5)
```

In plain words: the leaves (`LAI`) and the canopy efficiency (`ceff`) set how
much carbon the canopy could fix; this is then limited by how fast CO2 can enter
the leaves (`p_d`) and by how much light arrives (`E_0 I`), and scaled by how
long the day is (`h`).

Notes:
- `ceff` replaces the product "nitrogen-use efficiency x foliar nitrogen" of the
  original ACM, so foliar nitrogen is not needed.
- Temperature reaches GPP only through `e^(a8 Tmax)` (B3) and the range `Tr`
  (B1). Note B1 takes **half** the range; `Tr` itself is the full range.
- The `+ 10` in B5: Chuter et al. count days from 21 December, not 1 January.
  Without it, the whole season shifts 10 days late (M3).
- There is no water stress: `psi_d` and `R_tot` are fixed numbers.
- There is no cold limit on spring photosynthesis (Stettz et al. 2022,
  `DECISIONS.md` F9).

### 3.5 Not in the published model: frost cut-off

The previous repository added one rule: no photosynthesis on days when
`(Tmax + Tmin) / 2 < -2 degC`. Reason given: ACM was calibrated only between
about 7 and 30 degC, and Hyytiala reaches about -25 degC; Williams et al. (1997)
used a similar rule at their Oregon conifer sites. It is safe for the sampler
because it depends only on the inputs, never on a parameter. Whether to keep it
is open (M4).

## 4. Corrections to the published equations

Found by the previous repository by checking the paper against itself. To be
re-checked against the PDF.

1. **Allocation closure (Table 1, footnote 1).** The footnote gives
   `f_woo = 1 - f_auto - f_fol - f_lab`, leaving out `f_roo`. A1 to A4 and A10
   show five uses of GPP, so the correct closure is
   `f_auto + f_lab + f_fol + f_roo + f_woo = 1`. **Corrected.**
2. **A8 uses `cr_fall` where `d_fall` is expected.** A7 is anchored on
   `d_onset`; the matching form of A8 would be `doy - d_fall + psi_f`. As
   printed, `d_fall` has no effect on anything, so its posterior would just equal
   its prior. The previous repository kept the printed form. **Open (M5).**
3. **The text derivation writes `Phi_onset` where `Phi_fall` is meant.** Prose
   only; does not reach the code.
4. **The same sentence mixes `cr_onset` and `cr_fall`.** Prose only.

## 5. Outputs

Per day: GPP, Ra, Rh, Reco, NEE, LAI, and the six pools. Only NEE enters the
likelihood (on assimilable days). GPP and Reco from FLUXNET's partitioning are
used only for comparison, never fitted (`research_plan_v2.md`).

## 6. Parameters

Bloom and Williams (2015, Table 1) list 24 symbols. `f_woo` follows from the
closure, so there are **23 free parameters**: 17 process parameters and 6
starting pool sizes.

"B&W range" is the published prior range (a flat, "uniform" prior between the
two numbers). "Previous repository" says how the old code treated it.
**Scale** and **Prior** are filled in steps 2 and 3.

| # | Code name | Symbol | Meaning | Unit | B&W range | Previous repository | Scale | Prior |
|---|---|---|---|---|---|---|---|---|
| 1 | `f_auto` | f_auto | share of GPP breathed out by the trees | - | 0.3 to 0.7 | uniform, B&W range | | |
| 2 | `f_lab` | f_lab | share of GPP to labile | - | 0.01 to 0.5 | Dirichlet split of `1 - f_auto`, weights from site data | | |
| 3 | `f_fol` | f_fol | share of GPP to foliage | - | 0.01 to 0.5 | (Dirichlet, as above) | | |
| 4 | `f_roo` | f_roo | share of GPP to fine roots | - | 0.01 to 0.5 | (Dirichlet, as above) | | |
| - | `f_woo` | f_woo | share of GPP to wood | - | derived | (Dirichlet, as above) | derived | derived |
| 5 | `theta_roo` | theta_roo | fine-root turnover | d-1 | 1e-4 to 1e-2 | B&W range | | |
| 6 | `theta_woo` | theta_woo | wood turnover | d-1 | 2.5e-5 to 1e-3 | B&W range | | |
| 7 | `theta_lit` | theta_lit | litter respiration rate at 0 degC | d-1 | 1e-4 to 1e-2 | from a 1 to 5 year residence time | | |
| 8 | `theta_som` | theta_som | soil respiration rate at 0 degC | d-1 | 1e-7 to 1e-3 | from a 24 to 61 year residence time (site soil carbon) | | |
| 9 | `theta_min` | theta_min | litter-to-soil transfer at 0 degC | d-1 | 1e-5 to 1e-2 | B&W range | | |
| 10 | `temperature_exponent` | Theta | how fast decomposition rises with temperature | degC-1 | 0.018 to 0.08 | B&W range | | |
| 11 | `d_onset` | d_onset | day of leaf onset | day of year | 1 to 365 | uniform, B&W range (set by hand in one early run) | | |
| 12 | `cr_onset` | cr_onset | length of the leaf-onset pulse | days | 10 to 100 | uniform, B&W range | | |
| 13 | `d_fall` | d_fall | day of leaf fall | day of year | 1 to 365 | uniform, B&W range, but no effect as printed (M5) | | |
| 14 | `cr_fall` | cr_fall | length of the leaf-fall pulse | days | 20 to 150 | uniform, B&W range | | |
| 15 | `c_lf` | c_lf | share of needles shed per year | - | 1/8 to 1 | 0.20 to 0.33 (needles live 3 to 5 years) | | |
| 16 | `lma` | c_lma | leaf carbon per leaf area | g C m-2 | 10 to 400 | 116 to 247, from site litterfall and LAI | | |
| 17 | `ceff` | c_eff | canopy efficiency | - | 10 to 100 | 5 to 20 (REFLEX, Fox et al. 2009) | | |
| 18 | `c_lab_0` | C_lab(0) | starting labile pool | g C m-2 | 20 to 2000 | derived from site litterfall | | |
| 19 | `c_fol_0` | C_fol(0) | starting foliage pool | g C m-2 | 20 to 2000 | derived: litterfall / `c_lf` | | |
| 20 | `c_roo_0` | C_roo(0) | starting fine-root pool | g C m-2 | 20 to 2000 | derived from root litter and `theta_roo` | | |
| 21 | `c_woo_0` | C_woo(0) | starting wood pool | g C m-2 | 100 to 1e5 | derived: tree carbon minus the other plant pools | | |
| 22 | `c_lit_0` | C_lit(0) | starting litter pool | g C m-2 | 20 to 2000 | derived: two new parameters were sampled instead, total soil respiration at 0 degC (`rh_ref`) and the soil's share of it (`f_som`) | | |
| 23 | `c_som_0` | C_som(0) | starting soil pool | g C m-2 | 100 to 2e5 | derived (as above) | | |

**Why the previous repository moved away from the published ranges.** With the
published ranges, 837 of 1,000 draws from the prior produced physically
impossible fluxes at Hyytiala (for example, about 50 g C m-2 d-1 of soil
respiration, when the site's whole net exchange is around 0.6). The B&W ranges
were written to cover forests worldwide, so they are far too wide for one site.
This is exactly what the prior predictive check (Gelman 2.4) is meant to catch,
and it is why each prior needs its own decision here.

**Candidate scales (to decide in step 2).** Only a starting point:
- rates (`theta_...`) and pool sizes: log scale, because their ranges span 2 to
  4 powers of ten;
- allocation fractions: a simplex (a set of shares that always adds up to 1),
  which also fixes correction 1;
- days and periods: plain scale within their bounds.

## 7. Fixed ACM coefficients

From Chuter et al. (2015). They are fitted per site, not universal. The previous
repository used the **Loobos** set (evergreen forest, 52 N) as the closest
match; neither published set is boreal (M2). Values copied from the previous
repository and to be checked against Chuter et al., Appendix A.

| Coefficient | Loobos (used before) | Oregon pine (reference) | Role |
|---|---|---|---|
| a2 | 0.0156 | 0.0142 | day-length slope |
| a3 | 4.22273 | 0.980 | CO2 compensation point |
| a4 | 208.868 | 217.9 | CO2 half-saturation point |
| a5 | 0.0453 | 0.155 | day-length offset |
| a6 | 0.3783 | 2.653 | hydraulic resistance coefficient |
| a7 | 7.1929 | 4.309 | maximum light-use efficiency |
| a8 | 0.0111 | 0.060 | temperature coefficient (on Tmax) |
| a9 | 2.1001 | 1.062 | LAI coefficient of light-use efficiency |
| a10 | 0.7897 | 0.0006 | exponent on water potential |
| psi_d | 2.0 | 0.8502 | soil-to-leaf water potential difference (MPa) |
| R_tot | 1.0 | 1.0 | total hydraulic resistance |

## 8. Open model choices

To decide one by one before the priors. Each gets an entry in `DECISIONS.md`.

| ID | Question | Previous repository | Notes |
|---|---|---|---|
| M1 | Which ACM form: Chuter et al. (2015) or Williams et al. (1997)? | Chuter | The 1997 form has no latitude in its day length, which overstates winter GPP at 62 N by up to 12 times. Chuter's form is assumed to be what DALEC2 uses; the papers cannot fully confirm this. |
| M2 | Which ACM coefficient set? | Loobos | Neither set is boreal; a stated limitation either way. |
| M3 | Count days from 21 December in B5 (`doy + 10`)? | yes | Follows Chuter et al. Section II.A. |
| M4 | Keep the frost cut-off (no GPP below -2 degC mean)? | yes | Not in Bloom and Williams. Interacts with the Stettz et al. (2022) finding. |
| M5 | A8: keep `cr_fall` as printed, or correct to `d_fall`? | kept as printed | As printed, `d_fall` does nothing and must not be reported as "unidentifiable". |
| M6 | Phenology at an evergreen site: which of `d_onset`, `cr_onset`, `d_fall`, `cr_fall`, `c_lf` are fitted, and which fixed? | sampled over B&W ranges (`c_lf` narrowed to 3 to 5 year needle life) | NEE alone may say little about leaf timing in an evergreen stand (RQ1, RQ3). A range of 1 to 365 days for `d_onset` lets leaves appear in midwinter. |
| M7 | Starting pools: sample them, or derive them from site measurements? | derived (Ilvesniemi et al. 2009 site data) | Directly affects the prior predictive check. |
| M8 | Allocation: Dirichlet simplex? | yes, weights from site data | Listed as locked in the project notes; to be re-assessed like everything else. |
| M9 | LAI convention for `lma` (hemisurface or projected leaf area)? | both kept as a sensitivity run | Changes the `lma` range by about 25%. |

## 9. Still to check against the original papers

- Every equation in Section 3 against Bloom and Williams (2015) Appendix A and
  Chuter et al. (2015) Appendix B, including the A9 exponent (Section 3.3).
- The four corrections in Section 4.
- The ACM coefficients in Section 7 against Chuter et al. Appendix A and D.
- The previous repository's site numbers (litterfall, soil carbon, LAI) against
  Ilvesniemi et al. (2009) and Kolari (2010), before they are used in priors.

## References

- Bloom, A. A. and Williams, M. (2015). Constraining ecosystem carbon dynamics in a data-limited world. Biogeosciences 12, 1299-1315.
- Chuter, A. M., Aston, P. J., Skeldon, A. C. and Roulstone, I. (2015). A dynamical systems analysis of the data assimilation linked ecosystem carbon (DALEC) models. Chaos 25(3), 036401.
- Fox, A. et al. (2009). The REFLEX project. Agricultural and Forest Meteorology 149, 1597-1615.
- Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808.
- Ilvesniemi, H. et al. (2009). Long-term measurements of the carbon balance of a boreal Scots pine dominated forest ecosystem. Boreal Environment Research 14(4), 731-753. (Source of the previous repository's site carbon numbers, Fig. 6.)
- Kolari, P. (2010). Dissertationes Forestales 99, doi:10.14214/df.99. (Source of the previous repository's LAI numbers.)
- Stettz, S. G. et al. (2022). Resolving temperature limitation on spring productivity in an evergreen conifer forest using a model-data fusion framework. Biogeosciences 19, 541-558.
- Williams, M. et al. (1997). Predicting gross primary productivity in terrestrial ecosystems. Ecological Applications 7(3), 882-894.