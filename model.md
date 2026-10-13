# The DALEC2 model

What the model is, every equation it runs, and every parameter it has. It is the
starting point for the priors and the code.

This is **step 1** of Gelman et al. (2020), Section 2 ("Before fitting a model"):
write the model down before fitting it. Steps 2 and 3, choosing the scale and
prior for each parameter, fill the empty columns of the parameter table in
Section 6. The model-structure choices are decided (Section 8, `DECISIONS.md`
D21 to D30).

Status: 2026-10-13. Equations from Bloom and Williams (2015, Appendix A) and
Chuter et al. (2015, Appendix B), checked against the **official DALEC2 code**
(CARDAMOM, `DALEC_400`, the version used by Bloom et al. 2016) and the previous
repository. Where the printed paper and the code differ, this file follows the
code and says so. Equations must still be checked against the PDFs before they
go into the thesis (Section 9).

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
| `doy` | day of year (resets every 1 January) | day | from the date; used for day length |
| `t` | days since the start of the run (never resets) | day | from the date; used for leaf timing (D30) |

**Fixed site constant:** latitude 61.85 N (sets day length).

**Fixed model constants** (not fitted): the ACM coefficients `a2` to `a10`,
`psi_d` and `R_tot` (Section 7, D22), and the phenology constants `s`, 6.9088 and
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
temperature (D20); the official code computes it the same way.

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
`t` is the continuous day counter (D30), as in the official code.

```
(A7) Phi_onset = sqrt(2/pi) (6.9088 / cr_onset)
                 exp( -[ sin((t - d_onset - 0.6245 cr_onset) / s) sqrt(2) s / cr_onset ]^2 )

(A8) Phi_fall  = sqrt(2/pi) (-log(1 - c_lf) / cr_fall)
                 exp( -[ sin((t - d_fall + psi_f) / s) sqrt(2) s / cr_fall ]^2 )
                 (corrected: the paper prints cr_fall where d_fall belongs; D25)

     s = 365.25 / pi

(A9) psi solves  2 sqrt(pi) psi = log(1 - c_lf) exp(-psi^2),   psi_f = psi cr_fall / sqrt(2)
     (form recovered from the official code; computed with its polynomial, D31)
```

`psi` is not a parameter: it is calculated from `c_lf`. The official code uses
a fitted sixth-order polynomial (`offset.c`) in `log(L - 1)`, with leaf lifespan
`L = 1 / c_lf`:

```
psi = sum_k m_k (log(L - 1))^k,   k = 6 ... 0
m = 2.3599784710e-5, 3.32730053021e-4, 9.01865258885e-4, -5.437736864888e-3,
    -2.0836027517787e-2, 0.126972018064287, -0.188459767342504
```

**How A9 was settled (2026-10-13).** The polynomial was compared numerically
with candidate equations. It matches `2 sqrt(pi) psi = log(1 - c_lf) exp(-psi^2)`
to within 0.004 over leaf lifespans of 1.001 to 8 years. It does **not** match
the previous repository's form, `2 sqrt(pi) log(1 - c_lf) psi = exp(-psi^2)`,
which puts `log(1 - c_lf)` on the wrong side; at `c_lf = 0.25` that form gives
`psi = -0.646` against the code's `-0.079`, an offset about 8 times too large.
The `exp(-psi^2/2)` reading is not supported either. Two cross-checks: the same
polynomial at `L = 1.001` gives `-0.88318`, and `0.88318 / sqrt(2) = 0.6245`,
the constant in A7; and `log(1001) = 6.9088`, the other A7 constant. So A7 is
the onset version of the same rule. Whether the printed A9 has a typo or the
previous repository mis-copied it must be checked in the PDF.

The `sin` has a period of 730.5 days, and it is squared, so each pulse repeats
every 365.25 days.

Hyytiala is an **evergreen** pine forest: it keeps its needles all year and
sheds only the oldest ones. In DALEC2 terms, the onset pulse is the spring flush
of new shoots and the fall pulse is the autumn drop of old needles. All five
phenology parameters are fitted, with site-based ranges (D26).

### 3.4 Photosynthesis: the Aggregated Canopy Model (ACM)

Bloom and Williams (2015) give GPP only as a function from ACM (Williams et al.
1997). The form DALEC actually uses is written out by Chuter et al. (2015,
Appendix B, Eq. 12), which differs from the 1997 paper, most importantly in
using true day length from latitude. The official code (`ACM.c`) runs exactly
this form (D21).

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
  Without it, the whole season shifts 10 days late. The code does the same (D23).
- The code writes the declination as `-23.4 degrees x cos(...)`, which is the
  same as `-0.408 radians x cos(...)`.
- There is no water stress: `psi_d` and `R_tot` are fixed numbers.
- There is no cold limit on spring photosynthesis (Stettz et al. 2022,
  `DECISIONS.md` F9).

### 3.5 Optional switch, off by default: frost cut-off

Not part of DALEC2 (not in the paper, not in the official code). The previous
repository added it: no photosynthesis on days when `(Tmax + Tmin) / 2 < -2 degC`,
because ACM was calibrated only between about 7 and 30 degC and Hyytiala reaches
about -25 degC. The main model runs **without** it (D24). It stays available as
a switch for one prior predictive check and one sensitivity run. It is safe for
the sampler because it depends only on the inputs, never on a parameter.

## 4. Corrections to the published equations

Found by the previous repository by checking the paper against itself. To be
re-checked against the PDF.

1. **Allocation closure (Table 1, footnote 1).** The footnote gives
   `f_woo = 1 - f_auto - f_fol - f_lab`, leaving out `f_roo`. A1 to A4 and A10
   show five uses of GPP, so the closure used here is
   `f_auto + f_lab + f_fol + f_roo + f_woo = 1`. The official code splits the
   carbon step by step (foliage, then labile, then roots, wood takes the rest),
   so the footnote is best read as a shorthand of that, not a different model.
   **Handled by the Dirichlet split (D28).**
2. **A8 uses `cr_fall` where `d_fall` is expected.** A7 is anchored on
   `d_onset`; the official code anchors leaf fall on the leaf-fall day ("Fday").
   As printed, `d_fall` would have no effect. **Corrected (D25).**
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
**Scale** (step 2, D32) and **Prior** (step 3, D33 to D39) are now filled. All priors are smooth; `theta_min` is computed from the share `h` (D35).

| # | Code name | Symbol | Meaning | Unit | B&W range | Previous repository | Scale | Prior |
|---|---|---|---|---|---|---|---|---|
| 1 | `f_auto` | f_auto | share of GPP breathed out by the trees | - | 0.3 to 0.7 | uniform, B&W range | interval | Beta(25.5, 24.5), D34 |
| 2 | `f_lab` | f_lab | share of GPP to labile | - | 0.01 to 0.5 | Dirichlet split of `1 - f_auto`, weights from site data | simplex | (1 - f_auto) x Dirichlet, alpha 13.06, D34 |
| 3 | `f_fol` | f_fol | share of GPP to foliage | - | 0.01 to 0.5 | (Dirichlet, as above) | simplex | (Dirichlet) alpha 13.06, D34 |
| 4 | `f_roo` | f_roo | share of GPP to fine roots | - | 0.01 to 0.5 | (Dirichlet, as above) | simplex | (Dirichlet) alpha 15.26, D34 |
| - | `f_woo` | f_woo | share of GPP to wood | - | derived | (Dirichlet, as above) | simplex | (Dirichlet) alpha 44.27, D34 |
| 5 | `theta_roo` | theta_roo | fine-root turnover | d-1 | 1e-4 to 1e-2 | B&W range | log | LogNormal, median 1.78e-3, sdlog 0.45, D35 |
| 6 | `theta_woo` | theta_woo | wood turnover | d-1 | 2.5e-5 to 1e-3 | B&W range | log | LogNormal, median 5.33e-5, sdlog 0.39, D35 |
| 7 | `theta_lit` | theta_lit | litter respiration rate at 0 degC | d-1 | 1e-4 to 1e-2 | from a 1 to 5 year residence time | log | LogNormal, median 1.22e-3, sdlog 0.41, D35 |
| 8 | `theta_som` | theta_som | soil respiration rate at 0 degC | d-1 | 1e-7 to 1e-3 | from a 24 to 61 year residence time (site soil carbon) | log | LogNormal, median 7.15e-5, sdlog 0.24, D35 |
| 9 | `theta_min` | theta_min | litter-to-soil transfer at 0 degC | d-1 | 1e-5 to 1e-2 | B&W range | via h (interval) | derived: theta_lit h/(1-h), h ~ Beta(4, 16), D35 |
| 10 | `temperature_exponent` | Theta | how fast decomposition rises with temperature | degC-1 | 0.018 to 0.08 | B&W range | log | LogNormal, median 0.0336, sdlog 0.35, D36 |
| 11 | `d_onset` | d_onset | day of leaf onset | day of year | 1 to 365 | uniform, B&W range (set by hand in one early run) | plain | Normal(145, 15), D37 |
| 12 | `cr_onset` | cr_onset | length of the leaf-onset pulse | days | 10 to 100 | uniform, B&W range | log | LogNormal, median 40, sdlog 0.35, D37 |
| 13 | `d_fall` | d_fall | day of leaf fall | day of year | 1 to 365 | uniform, B&W range, but had no effect (printed A8); active now (D25) | plain | Normal(265, 15), D37 |
| 14 | `cr_fall` | cr_fall | length of the leaf-fall pulse | days | 20 to 150 | uniform, B&W range | log | LogNormal, median 60, sdlog 0.35, D37 |
| 15 | `c_lf` | c_lf | share of needles shed per year | - | 1/8 to 1 | 0.20 to 0.33 (needles live 3 to 5 years) | interval | Beta(43.84, 122.84), D43 (round 3; was Beta(14.51, 39.24), D37) |
| 16 | `lma` | c_lma | leaf carbon per leaf area | g C m-2 | 10 to 400 | 116 to 247, from site litterfall and LAI | log | LogNormal, median 185, sdlog 0.25, D38 |
| 17 | `ceff` | c_eff | canopy efficiency | - | 10 to 100 | 5 to 20 (REFLEX, Fox et al. 2009) | log | LogNormal, median 10, sdlog 0.354, D42 (round 2; was median 16.5, D38) |
| 18 | `c_lab_0` | C_lab(0) | starting labile pool | g C m-2 | 20 to 2000 | derived from site litterfall | log | LogNormal, median 77, sdlog 0.50, D39 |
| 19 | `c_fol_0` | C_fol(0) | starting foliage pool | g C m-2 | 20 to 2000 | derived: litterfall / `c_lf` | log | LogNormal, median 590, sdlog 0.30, D39 |
| 20 | `c_roo_0` | C_roo(0) | starting fine-root pool | g C m-2 | 20 to 2000 | derived from root litter and `theta_roo` | log | LogNormal, median 238, sdlog 0.30, D39 |
| 21 | `c_woo_0` | C_woo(0) | starting wood pool | g C m-2 | 100 to 1e5 | derived: tree carbon minus the other plant pools | log | LogNormal, median 5900, sdlog 0.25, D39 |
| 22 | `c_lit_0` | C_lit(0) | starting litter pool | g C m-2 | 20 to 2000 | derived: two new parameters were sampled instead, total soil respiration at 0 degC (`rh_ref`) and the soil's share of it (`f_som`) | log | LogNormal, median 450, sdlog 0.50, D39 |
| 23 | `c_som_0` | C_som(0) | starting soil pool | g C m-2 | 100 to 2e5 | derived (as above) | log | LogNormal, median 6100, sdlog 0.25, D39 |

**Why the previous repository moved away from the published ranges.** With the
published ranges, 837 of 1,000 draws from the prior produced physically
impossible fluxes at Hyytiala (for example, about 50 g C m-2 d-1 of soil
respiration, when the site's whole net exchange is around 0.6). The B&W ranges
were written to cover forests worldwide, so they are far too wide for one site.
This is exactly what the prior predictive check (Gelman 2.4) is meant to catch,
and it is why each prior needs its own decision here.

**The official code uses some different ranges** (Bloom et al. 2016 version):
`f_auto` 0.2 to 0.8; `f_fol` 0.01 to 0.5 and `f_lab` 0.01 to 0.5 as shares of
NPP, `f_roo` 0.01 to 1 of what remains; leaf lifespan 1.001 to 8 years; `ceff`
5 to 50; `lma` 5 to 200 g C m-2 (from Kattge et al. 2011); onset and fall days
365.25 to 1461 on the continuous counter; starting pools from 1 g C m-2. Useful
in step 3.

**Scales (step 2, decided 2026-10-13, D32).** The paper gives each parameter a
*range*; the *scale* is how the sampler moves inside that range. It is read off
the range with one rule:

- **log**: the range covers a factor of 10 or more (upper / lower >= 10). The
  sampler works on `log(value)`, so 0.0001 and 0.01 are as easy to move between
  as 1 and 100.
- **interval**: the range covers less than a factor of 10, or the parameter is a
  position (a day of the year). The sampler works on an unbounded stretched
  version of the range (the "logit"), so it can never step outside it.
- **simplex**: the allocation split (D28), whose shares always add up to 1.

The rule is applied again to the final step-3 ranges: if a range is narrowed
below a factor of 10 (for example `lma`), its scale may change to interval.

## 7. Fixed ACM coefficients (D22)

Fixed inside DALEC2 itself: they are the same numbers in every run of the
official code (`DALEC_400.c`). They equal the "Loobos" set of Chuter et al.
(2015, Appendix A; evergreen forest, 52 N). The previous repository used the
same set, rounded. They were not fitted for a boreal forest: stated limitation.

| Coefficient | DALEC2 value (used) | Previous repository (rounded) | Role |
|---|---|---|---|
| a2 | 0.0156935 | 0.0156 | day-length slope |
| a3 | 4.22273 | 4.22273 | CO2 compensation point |
| a4 | 208.868 | 208.868 | CO2 half-saturation point |
| a5 | 0.0453194 | 0.0453 | day-length offset |
| a6 | 0.37836 | 0.3783 | hydraulic resistance coefficient |
| a7 | 7.19298 | 7.1929 | maximum light-use efficiency |
| a8 | 0.011136 | 0.0111 | temperature coefficient (on Tmax) |
| a9 | 2.1001 | 2.1001 | LAI coefficient of light-use efficiency |
| a10 | 0.789798 | 0.7897 | exponent on water potential |
| psi_d | -2.0 (MPa) | 2.0 | soil-to-leaf water potential difference; B1 uses its size |
| R_tot | 1.0 | 1.0 | total hydraulic resistance |

## 8. Model-structure choices (decided 2026-10-13)

| ID | Question | Decision | Entry |
|---|---|---|---|
| M1 | ACM form | Chuter et al. (2015) form, as in the official code | D21 |
| M2 | ACM coefficients | DALEC2's fixed constants, full precision | D22 |
| M3 | Day count for the declination | from 21 December (`doy + 10`) | D23 |
| M4 | Frost cut-off | off in the main model; switch for one sensitivity run | D24 |
| M5 | A8 anchor | corrected to `d_fall` | D25 |
| M6 | Phenology parameters | all five fitted, site-based ranges (step 3) | D26 |
| M7 | Starting pools | sampled, priors centred on site measurements (step 3) | D27 |
| M8 | Allocation | `f_auto` plus a Dirichlet split | D28 |
| M9 | LAI convention for `lma` | decided with the `lma` prior (step 3) | D29 |
| M10 | Time counter for phenology | continuous day counter `t` | D30 |

## 9. Still to check against the original papers

- Every equation in Section 3 against Bloom and Williams (2015) Appendix A and
  Chuter et al. (2015) Appendix B, including the A9 exponent (Section 3.3).
  The structure already matches the official code; the check is for the
  thesis text and citations.
- The four corrections in Section 4.
- The ACM coefficients in Section 7 against Chuter et al. Appendix A (they
  already match the code).
- The previous repository's site numbers (litterfall, soil carbon, LAI) against
  Ilvesniemi et al. (2009) and Kolari (2010), before they are used in priors.

## References

- Bloom, A. A., Exbrayat, J.-F., van der Velde, I. R., Feng, L. and Williams, M. (2016). The decadal state of the terrestrial carbon cycle: global retrievals of terrestrial carbon allocation, pools, and residence times. PNAS 113(5), 1285-1290.
- Bloom, A. A. and Williams, M. (2015). Constraining ecosystem carbon dynamics in a data-limited world. Biogeosciences 12, 1299-1315.
- CARDAMOM framework, DALEC2 source code: github.com/CARDAMOM-framework/CARDAMOM, `C/projects/CARDAMOM_MODELS/DALEC/DALEC_400/` and `DALEC_ALL/ACM.c`, `offset.c` (read 2026-10-13).
- Chuter, A. M., Aston, P. J., Skeldon, A. C. and Roulstone, I. (2015). A dynamical systems analysis of the data assimilation linked ecosystem carbon (DALEC) models. Chaos 25(3), 036401.
- Fox, A. et al. (2009). The REFLEX project. Agricultural and Forest Meteorology 149, 1597-1615.
- Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808.
- Ilvesniemi, H. et al. (2009). Long-term measurements of the carbon balance of a boreal Scots pine dominated forest ecosystem. Boreal Environment Research 14(4), 731-753. (Source of the previous repository's site carbon numbers, Fig. 6.)
- Kolari, P. (2010). Dissertationes Forestales 99, doi:10.14214/df.99. (Source of the previous repository's LAI numbers.)
- Stettz, S. G. et al. (2022). Resolving temperature limitation on spring productivity in an evergreen conifer forest using a model-data fusion framework. Biogeosciences 19, 541-558.
- Williams, M. et al. (1997). Predicting gross primary productivity in terrestrial ecosystems. Ecological Applications 7(3), 882-894.