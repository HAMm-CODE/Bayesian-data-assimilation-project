MASTER'S THESIS: RESEARCH PLAN (version 3)
Bayesian Data Assimilation of Forest Ecosystem Carbon Dynamics Using Eddy-Covariance Data
Tampere University
Author: Hamza Sentongo
Updated: 1 November 2026 (replaces version 2 of 4 July 2026)

Note on this version: it records what was actually decided and done, so that it can be used as the reference for writing. Every choice points to its entry in DECISIONS.md (D = decision, F = finding, L = likelihood). It is a planning document, not thesis text.

---

1. Aim

The thesis asks how far eddy-covariance observations of net ecosystem exchange (NEE) constrain the parameters and the carbon dynamics of a process-based forest carbon model, DALEC2, at one boreal site (Hyytiälä, FI-Hyy), within a Bayesian framework. The analysis follows the Bayesian workflow of Gelman et al. (2020).

2. Research questions

Main question: To what extent do eddy-covariance observations of NEE constrain the parameters and carbon dynamics of a process-based forest carbon model within a Bayesian framework?

- RQ1. To what extent can NEE observations from a single boreal site constrain the latent carbon pools and the key physiological parameters of the model?
- RQ2. How does uncertainty propagate from the observations to the model parameters, pools and fluxes?
- RQ3. Which parameters are identifiable from flux data alone, and which remain poorly constrained or confounded with one another?

To do: settle the final wording of RQ3 (an earlier version also mentioned what gradient-based diagnostics reveal). The measures used for each RQ are listed in Section 7.

3. Site and data (D4, D5, D12, D15, D18 to D20, F3 to F6)

- Site: Hyytiälä SMEAR II, FI-Hyy (61.85 N, 24.29 E, about 181 m), a managed Scots pine stand sown in 1962 (Hari and Kulmala 2005).
- Data: FLUXNET2015 FULLSET, daily file (NEE, its quality flag and its random uncertainty RANDUNC, shortwave radiation, CO2) and half-hourly file (air temperature, used for the daily maximum and minimum) (Pastorello et al. 2020).
- Model inputs (D19, D20): daily maximum and minimum air temperature, shortwave radiation (MJ m-2 d-1), atmospheric CO2, day of year, and the site latitude. No precipitation or VPD: DALEC2 has no water cycle.
- Years: fitting 1997 to 2010; testing 2011 to 2014 (D4).
- Days used in the likelihood (D12): NEE present, quality flag at least 0.75, and RANDUNC present and above zero. This gives 4,413 of 5,113 days in 1997 to 2010 and 810 test days. RANDUNC is missing for all of 2014 (F5), so in practice the test days are in 2011 to 2013. No day is deleted; the model runs through every day.
- GPP and ecosystem respiration from FLUXNET's partitioning are used only for comparison, never fitted.

4. The model (D20 to D31, D44, D55, F33)

- DALEC2 (Bloom and Williams 2015), part of the DALEC family (Williams et al. 2005): six carbon pools (labile, foliage, fine roots, wood, litter, soil organic matter), daily time step, 17 process parameters and 6 starting pool sizes.
- Photosynthesis: the Aggregated Canopy Model (ACM; Williams et al. 1997). Its coefficients are fixed at the published DALEC2 values (not fitted for a boreal site); only canopy efficiency (c_eff) is estimated. LAI = foliage carbon / LMA, and LAI feeds back into GPP every day.
- Implementation follows the official DALEC2 code (CARDAMOM `DALEC_400`; Bloom et al. 2016), confirmed line by line by a second official release (Williams 2019). Where the printed equations differ from the code (leaf fall anchored on the leaf-fall day, A8), the code is followed (F33).
- Extension: a cold limit on photosynthesis (Stettz et al. 2022; code `DALEC_831`): GPP is multiplied by a ramp on the daily minimum temperature, 0 at or below T_shutdown and 1 at or above T_full. Two extra parameters, so 25 in total.
- Ecological constraints: one. A smooth version of the steady-state constraint of the official code (EDC2_400.c; Bloom and Williams 2015), applied to the foliage pool only (D44). It lowers the prior probability of parameter sets whose needle pool runs away or collapses over the 14 years. The official code applies a hard version to all pools; the smooth form is used because NUTS cannot work with hard limits, and foliage only because the stand is growing and was thinned.
- Two programs give the same numbers: a plain numpy version (readable reference) and a PyTensor version (for NUTS gradients), checked against each other by tests (D46).

5. Bayesian set-up

5.1 Priors (D32 to D39, D42, D43, D55, D66, D67)
- 25 sampled quantities, the same as the model's 25 parameters.
- Every prior is smooth and centred on a Hyytiälä measurement or a literature value (mainly Ilvesniemi et al. 2009 and Kolari 2010), wide enough to keep the published range possible (D33). The official DALEC2 priors are flat between limits; smooth shapes are used because hard limits are walls for NUTS.
- Shapes: log-normal for rates, pools, LMA and c_eff; beta for fractions; normal for the leaf-onset and leaf-fall days.
- Allocation (D67): f_auto (share of GPP respired) has a beta prior; the rest is split step by step, as in the official code: foliage takes a fraction of the growth carbon, labile a fraction of what is left, fine roots a fraction of what is left after that, wood the rest; each fraction has a beta prior.
- Litter to soil (D66): the rate theta_min gets a log-normal prior whose 95% range is the published range of Bloom and Williams (2015), 1e-5 to 1e-2 per day.
- Cold thresholds (D55): each a Beta(2, 2) fraction of the official ranges (-15 to 0 degC and 0 to 15 degC).
- Two priors were revised after the prior predictive check, each with a stated reason: c_eff (D42) and needle life c_lf (D43).

5.2 Likelihood (L1, L2)
- NEE_t ~ Normal(model NEE_t, sqrt(RANDUNC_t^2 + sigma_m^2)) on the assimilable days, with sigma_m ~ HalfNormal(1).
- sigma_m is the model error: RANDUNC covers measurement error only. Planned as a sensitivity option before the data were used; adopted after the first real-data fits, where RANDUNC alone gave separate modes and impossible forests (F17, F18).
- Assumptions stated as limitations: Gaussian errors (Hollinger and Richardson 2005 suggest heavier tails); independent days, although neighbouring days are correlated (F11); one sigma_m for the whole year.

6. Computation and workflow

6.1 Sampling (D46, D47, D49, D50, D53)
- NUTS in PyMC, with the model in PyTensor for automatic gradients; 4 chains of 1,000 warm-up and 1,000 kept draws; dense mass matrix; target_accept 0.99; each chain starts from a physically possible prior draw. Runs on CSC Roihu (project_2020170).
- Computation checks: no divergences, R-hat below 1.01, effective sample size at least 400 (Vehtari et al. 2021).

6.2 Workflow steps, each with a pass rule fixed before it was run (Gelman et al. 2020)
1. Prior predictive check (D40, D57): 1,000 prior draws run with the real weather, no NEE. Tier 1 (impossible forests) at most 5% (7% for the cold model); Tier 2 (plausible for the site) at least 50%. The soft foliage constraint is included by weighting the draws.
2. Fake-data check (D48, D49, D63): a "true" forest (the first fully plausible prior draw) produces fake NEE with real-size noise on the real days; the fit must pass the computation checks and contain the truth in its 90% interval for all but a few quantities (the number allowed follows from chance: at most 12% of perfect fits would miss more).
3. Real-data fit and posterior checks (D51): computation checks first; then prior against posterior width for every parameter; the forest the posterior implies against site values; posterior predictive coverage overall and by season.
4. Held-back years (D61): 1,000 posterior draws run from 1997 through 2014. Decision score: log predictive density on the test days; "clearly better" means a difference above twice its standard error, widened for correlated days. PSIS-LOO on the fitted years is reported only.

7. Model development and model choice (D55, D62 to D65, F27 to F36)

- The original DALEC2 (no cold limit) could not be fitted to a believable forest on real NEE: separate modes, no believable chain (F27, F30). It is reported, not fitted again (D62). The cold-limit model is the base model (M2).
- Supervisor's request: a simpler model, holding constant what the fake data show NEE cannot learn. Ladder of reduced models (D63), chosen only from the fake-data check: M2 (25 free), MA (starting wood pool fixed), MB (+ starting soil pool), MC (+ soil and wood turnover), M0 (everything not strongly informed, 14 free).
- Choice rule (D63): a model takes part only if its fake-data check and its fit pass; among models within twice the widened standard error of the best test score, the smallest wins; a believable forest is preferred.
- Literature posteriors as priors (D64, D65, F35): the only independent DALEC2 posterior for this location (CARDAMOM, Bloom et al. 2016) describes a 1 degree, deciduous-like, peat-soil ecosystem (LMA 23, needle life 0.4 years, soil carbon four times the stand). It is not used as a prior for the stand; it is compared in the discussion. It agrees with this thesis in one striking way: a thin canopy and high heterotrophic respiration.
- Result with the earlier priors (F34): MA chosen; no model gives a believable forest (peak LAI about 1, Rh about 555 against the site's 290 to 370).
- Current status (1 November 2026, F36): with the new priors (D66, D67), no model is yet eligible. The fake-data misses come mainly from an unusual "true" forest in the tails of the priors; theta_min mixes slowly. Proposed next step (D68, to be confirmed): a second true forest with the coverage judged over both truths, more draws and slightly smaller steps, M0 not rerun.

8. How each research question is answered

- RQ1: posterior width divided by prior width for each parameter (below 0.5 strongly informed, 0.5 to 0.9 partly, 0.9 or more not informed), on fake data and on real data.
- RQ2: spread of the pools and fluxes across posterior draws; the parts of the predictive spread due to parameters, measurement error and model error; coverage on the test years.
- RQ3: posterior correlations between parameters; the reduced-model ladder (what can be fixed without loss); the thin-canopy trade-off between canopy size, leaf efficiency and allocation, and its link to the fixed ACM coefficients (F22, F28); the CARDAMOM comparison (F35).

9. Expected results (updated)

(i) Posterior distributions for DALEC2 with a cold limit at FI-Hyy from NEE alone, with a clear split between strongly, partly and not informed parameters.
(ii) Quantified uncertainty for the fluxes and pools, and its test-year performance.
(iii) Evidence that NEE alone can produce a model that predicts NEE reasonably well from a forest whose inside is wrong (thin canopy, high Rh), and a discussion of why.

10. Preliminary table of contents

1 Introduction (background and motivation; research questions; scope and structure)
2 Theoretical background (forest carbon and NEE; eddy covariance and FLUXNET2015; the DALEC family; Bayesian data assimilation and MCMC; identifiability and equifinality)
3 Materials and methods
  3.1 Site and data
  3.2 The forward model
  3.3 Bayesian set-up (3.3.1 priors, 3.3.2 likelihood)
  3.4 Workflow (prior check, fake-data check, real fit and checks, held-back years)
  3.5 Sampling and computation
  3.6 Model development and model choice
  3.7 Software and reproducibility
4 Results (prior and fake-data checks; real-data fits; RQ1; RQ2; RQ3; model choice)
5 Discussion (interpretation and comparison with earlier studies; limitations; further work)
6 Conclusions
Appendix A: model equations, and the differences between the printed equations and the official code
Use of AI page, after the abstract

11. Still to do

- Settle D68 and rerun; score the test years; apply the choice rule.
- Figures and tables for the chosen model (RQ1 to RQ3).
- Freeze the software versions (D2).
- Small checks (OPEN-3: one definition of spring onset; OPEN-6: references marked "check", and tell the supervisor about the soft foliage constraint).
- Writing: chapters 1 to 6, appendix, abstract, AI-use page and log.

References (as used so far; check each against the original before citing)

Bloom, A. A. and Williams, M. (2015). Constraining ecosystem carbon dynamics in a data-limited world: integrating ecological "common sense" in a model-data fusion framework. Biogeosciences 12, 1299-1315.
Bloom, A. A., Exbrayat, J.-F., van der Velde, I. R., Feng, L. and Williams, M. (2016). The decadal state of the terrestrial carbon cycle: global retrievals of terrestrial carbon allocation, pools, and residence times. PNAS 113(5), 1285-1290.
CARDAMOM framework, source code. github.com/CARDAMOM-framework/CARDAMOM (DALEC_400, DALEC_831).
Famiglietti, C. A. et al. (2021). Optimal model complexity for terrestrial carbon cycle prediction. Biogeosciences 18, 2727.
Fox, A. et al. (2009). The REFLEX project: comparing different algorithms and implementations for the inversion of a terrestrial ecosystem model against eddy covariance data. Agricultural and Forest Meteorology 149, 1597-1615.
Gelman, A. et al. (2020). Bayesian Workflow. arXiv:2011.01808.
Hari, P. and Kulmala, M. (2005). Station for Measuring Ecosystem-Atmosphere Relations (SMEAR II). Boreal Environment Research 10(5), 315-322.
Hollinger, D. Y. and Richardson, A. D. (2005). Uncertainty in eddy covariance measurements and its application to physiological models. Tree Physiology 25.
Ilvesniemi, H. et al. (2009). Long-term measurements of the carbon balance of a boreal Scots pine dominated forest ecosystem. Boreal Environment Research 14(4), 731-753.
Kolari, P. (2010). Dissertationes Forestales 99, doi:10.14214/df.99.
Pastorello, G. et al. (2020). The FLUXNET2015 dataset and the ONEFlux processing pipeline for eddy covariance data. Scientific Data 7, 225.
Stettz, S. G. et al. (2022). Resolving temperature limitation on spring productivity in an evergreen conifer forest using a model-data fusion framework. Biogeosciences 19, 541-558.
Vehtari, A., Gelman, A. and Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing 27(5), 1413-1432.
Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C. (2021). Rank-normalization, folding, and localization: an improved R-hat for assessing convergence of MCMC. Bayesian Analysis 16(2), 667-718.
Williams, M. et al. (1997). Predicting gross primary productivity in terrestrial ecosystems. Ecological Applications 7(3), 882-894.
Williams, M., Schwarz, P. A., Law, B. E., Irvine, J. and Kurpius, M. R. (2005). An improved analysis of forest carbon dynamics using data assimilation. Global Change Biology 11(1), 89-105.
Williams, M. (2019). DALEC2 [software]. University of Edinburgh. https://doi.org/10.7488/ds/2504.
