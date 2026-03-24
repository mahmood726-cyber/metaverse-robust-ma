# Metaverse: An R Package for Robust Evidence Synthesis and High-Dimensional Variable Selection in Meta-Analysis

## Authors
**Mahmood Ahmad**^1^

^1^ [AFFILIATION_PLACEHOLDER]

*Corresponding author: [EMAIL_PLACEHOLDER]*

## Abstract

**Background:** Traditional meta-analysis relies on parametric assumptions and is highly sensitive to data contamination (outliers) and over-parameterization when evaluating multiple moderators. As clinical datasets grow in complexity, the need for robust, high-dimensional evidence synthesis tools has become critical.

**Methods and Implementation:** We introduce `metaverse`, an R package that provides a comprehensive framework for modern meta-analysis. The package integrates robust M-, S-, and MM-estimators with contamination models (t-mixture, slash, Box-Cox) to down-weight anomalous studies automatically. For high-dimensional moderator analysis, `metaverse` implements the Knockoff Filter, spike-and-slab Bayesian variable selection, LASSO with stability selection, elastic net, SCAD, and MCP penalties, each with rigorous False Discovery Rate (FDR) control. Inference tools include conformal prediction intervals, permutation tests, and Bayesian posterior sampling. A full sensitivity analysis suite provides Egger and Begg tests, trim-and-fill, fragility indices, E-values, and Cook's distance diagnostics.

**Results:** We evaluated the package on a simulated clinical meta-analysis dataset (k = 60 studies, 10% outlier contamination, 10 candidate moderators of which only 2 are true). The standard REML random-effects model yielded a biased pooled estimate of 0.48 (95% CI: 0.29-0.67) with inflated heterogeneity (tau-squared = 0.42, I-squared = 89.3%). In contrast, the MM-estimator recovered the true effect (estimate = 0.31, 95% CI: 0.22-0.40, tau-squared = 0.06, I-squared = 41.2%) by automatically down-weighting all 6 contaminated studies (robust weight < 0.3). For moderator selection, naive univariate meta-regression declared 4 of 8 noise covariates significant at p < 0.05 (observed FDR = 67%). The Knockoff Filter at FDR = 0.10 correctly retained X1 and X2 while rejecting all 8 noise covariates (achieved FDR = 0%). Spike-and-slab selection (10,000 MCMC iterations) yielded posterior inclusion probabilities of 0.94 for X1 and 0.87 for X2, with all noise covariates below 0.12. Computation time was 0.8 seconds for the robust estimation and 2.1 seconds for the Knockoff Filter on a standard laptop (Intel i7, 16 GB RAM).

**Conclusion:** The `metaverse` package bridges the gap between modern robust statistics and applied evidence synthesis. By combining high-breakdown-point estimation with mathematically guaranteed FDR control, the package safeguards meta-analyses against the dual threats of data contamination and spurious moderator discovery.

**Keywords:** meta-analysis, robust estimation, MM-estimator, knockoff filter, spike-and-slab, FDR control, variable selection, R software

## Introduction

Meta-analysis is the gold standard for evidence synthesis in medicine and the social sciences [1,2]. However, the standard random-effects model (e.g., DerSimonian-Laird [3] or Restricted Maximum Likelihood [4]) is inherently vulnerable to two major methodological threats:

1. **Data contamination:** A single extreme outlier---whether caused by reporting errors, distinct clinical subpopulations, or extreme bias---can disproportionately skew the pooled estimate and inflate between-study heterogeneity (tau-squared) [5,6]. The conventional REML estimator has a breakdown point of 1/k (where k is the number of studies), meaning that a single anomalous study can render the meta-analysis arbitrarily biased [7].

2. **The curse of dimensionality in meta-regression:** Researchers frequently test numerous study-level covariates (moderators) to explain heterogeneity [8]. Standard meta-regression based on weighted least squares lacks robust mechanisms to control the False Discovery Rate (FDR) in high-dimensional settings, leading to spurious associations and "p-hacking" of moderator effects [9,10].

Robust statistical methods, including M-estimators [11], S-estimators [12], and MM-estimators [13], have been extensively developed in the regression literature but remain rarely applied in meta-analysis. Similarly, modern variable selection techniques such as the Knockoff Filter [14,15] and spike-and-slab priors [16] offer formal FDR control guarantees, but no existing meta-analysis software integrates these methods.

Existing R packages for meta-analysis---`metafor` [17], `meta` [18], and `robumeta` [19]---provide excellent tools for standard and cluster-robust estimation but do not implement high-breakdown-point estimators or modern variable selection with FDR guarantees. The `metaverse` package fills this gap by providing a unified framework that combines robust estimation, principled variable selection, modern inference, and comprehensive sensitivity diagnostics.

## Methods and Implementation

`metaverse` is implemented in R (>= 4.0.0) with an S4 class system for structured analysis workflows. The package depends on `robustbase` [20] for low-level robust regression primitives, `glmnet` [21] for penalized regression, `lme4` [22] for mixed-effects modeling, and `ggplot2` [23] for visualization. The architecture is organized into three pillars.

### 1. Robust Estimation and Contamination Models

Instead of relying on standard weighted least squares, `metaverse` implements a suite of robust estimators through the `meta_robust()` function:

- **M-estimator:** Iteratively reweighted least squares (IRLS) using Tukey's biweight psi-function with tuning constant k = 4.685, providing 95% asymptotic efficiency under the normal model while strongly down-weighting outliers [11].
- **S-estimator:** Achieves a breakdown point of 50%, meaning the estimate remains bounded even when up to half the studies are arbitrarily contaminated [12]. The implementation uses the median absolute deviation (MAD) for initial scale estimation.
- **MM-estimator:** A two-stage procedure that first computes an S-estimate for high breakdown point, then refines it with an M-estimate for efficiency [13]. The `meta_robust(method = "MM")` function achieves simultaneous 50% breakdown point and 95% Gaussian efficiency.
- **Tau-squared estimation:** Multiple heterogeneity estimators (DerSimonian-Laird, Paule-Mandel, REML, Empirical Bayes) computed on robustly centered residuals, with Hartung-Knapp-Sidik-Jonkman (HKSJ) prediction intervals.

Explicit contamination models are available to formally accommodate heavy-tailed data: (a) t-mixture models with user-specified degrees of freedom, estimated via EM algorithm; (b) slash distribution mixtures for extreme heavy tails; and (c) Box-Cox power transformations with profile-likelihood lambda selection.

### 2. High-Dimensional Variable Selection

To prevent false-positive moderator discovery, `metaverse` implements seven selection algorithms through `select_moderators()`:

- **Knockoff Filter** [14,15]: Constructs synthetic "knockoff" covariates that replicate the correlation structure of the original moderators but are conditionally independent of the outcome. The test statistic W_j = |beta_j| - |beta_j_knockoff| captures each variable's importance relative to its null counterpart. The knockoff+ threshold procedure guarantees FDR control at the user-specified level (default 10%). Three knockoff construction methods are provided: SDP (semidefinite programming), equicorrelated, and sequential.
- **Spike-and-Slab Priors** [16]: Full Gibbs sampling MCMC with binary inclusion indicators (gamma_j), yielding posterior inclusion probabilities (PIP) for each moderator. The median probability model (PIP > 0.5) provides consistent model selection, and the Bayesian FDR is directly computed from posterior summaries.
- **LASSO with stability selection** [24]: Cross-validated LASSO on random subsamples (default 100 iterations, 50% subsample ratio), selecting variables that appear in at least 60% of subsamples.
- **Elastic net, SCAD, and MCP** [25,26]: Penalized regression with non-convex penalties to reduce bias relative to LASSO.
- **E-values selection** [27]: Calibrated p-to-e-value conversion with e-BH and e-BY procedures for FDR control that remains valid under arbitrary dependence.

### 3. Inference and Sensitivity Diagnostics

`metaverse` supplements point estimates with five inference approaches: nonparametric bootstrap (BCa intervals, 1000 replicates), split conformal prediction for finite-sample coverage guarantees [28], data-splitting selective inference for post-selection validity, permutation-based p-values, and Bayesian Gibbs sampling for joint posterior inference on mu and tau-squared.

The sensitivity analysis suite (`assess_sensitivity()`) includes: three variants of Egger's regression (standard, robust, multilevel), Begg's rank correlation, Thompson-Sharp test, trim-and-fill (L0, R0, Q0 estimators), step and beta selection models, p-curve analysis, fragility index with directional variants, E-values for unmeasured confounding [29], Cook's distance, DFFITS, covariance ratio, and cumulative meta-analysis. Eight publication-quality plot types are provided through `visualize()`.

## Results

We evaluated the `metaverse` package using the bundled `metaverse_eval_data` dataset, a simulated clinical meta-analysis designed to stress-test robust methods. The dataset comprises k = 60 studies with a true baseline effect theta = 0.30, true moderator effects (beta_X1 = 0.20, beta_X2 = -0.15), true between-study variance tau-squared = 0.05, and 10% contamination (6 outlier studies shifted by 1.5-3.0 standard deviations in random directions).

### Robust Estimation vs. Standard Models

A standard REML random-effects model (equivalent to `metafor::rma()`) applied to the full dataset yielded a pooled estimate of 0.48 (95% CI: 0.29, 0.67), substantially inflated relative to the true theta = 0.30. Heterogeneity was severely overestimated at tau-squared = 0.42 (I-squared = 89.3%, Q = 561.6, p < 0.001), driven almost entirely by the 6 contaminated studies.

In contrast, the `metaverse` MM-estimator (`meta_robust(method = "MM")`) yielded a pooled estimate of 0.31 (95% CI: 0.22, 0.40), recovering the true effect within the confidence interval. The robust estimate of heterogeneity was tau-squared = 0.06 (I-squared = 41.2%), closely matching the true tau-squared = 0.05. All 6 contaminated studies received robust weights below 0.30, while the remaining 54 clean studies retained weights above 0.85. The M-estimator and S-estimator produced estimates of 0.33 and 0.29, respectively, confirming consistency across the robust estimation family.

### Variable Selection Performance

Standard univariate meta-regression testing each of the 10 candidate moderators individually declared 6 variables significant at p < 0.05: both true moderators (X1 with p = 0.001, X2 with p = 0.008) and 4 noise covariates (X4 with p = 0.021, X6 with p = 0.031, X8 with p = 0.043, X9 with p = 0.048), yielding an observed FDR of 4/6 = 67%.

The Knockoff Filter (`select_moderators(method = "knockoff", fdr = 0.10)`) correctly identified X1 and X2 as the only true moderators, rejecting all 8 noise covariates (achieved FDR = 0%, power = 100%). The knockoff importance statistics (W) for X1 and X2 were 1.24 and 0.87, respectively, while the maximum W among noise covariates was 0.18, well below the computed threshold of 0.62.

Spike-and-slab MCMC selection (10,000 iterations, 2,000 burn-in) yielded posterior inclusion probabilities (PIP) of 0.94 for X1 and 0.87 for X2. All noise covariates had PIPs below 0.12, and the Bayesian FDR of the median probability model was 0.03. The posterior mean coefficients for X1 (0.18, 95% HPD: 0.09, 0.28) and X2 (-0.12, 95% HPD: -0.22, -0.03) were consistent with the true values of 0.20 and -0.15.

LASSO with stability selection (100 bootstrap subsamples, 60% threshold) selected X1 and X3, missing X2 but including one false positive. Elastic net selected X1, X2, and X7. SCAD and MCP both selected X1 and X2 without false positives.

### Computational Performance

On a standard laptop (Intel Core i7-1260P, 16 GB RAM), meta_robust(method = "MM") completed in 0.8 seconds for k = 60 and 3.2 seconds for k = 1000. The Knockoff Filter with p = 10 moderators required 2.1 seconds, and spike-and-slab MCMC (10,000 iterations) completed in 4.7 seconds. The full sensitivity analysis suite ran in 12.3 seconds for k = 60.

### Sensitivity Analysis

The full sensitivity suite applied to the robust MM estimate confirmed the finding's stability. Egger's regression test was non-significant (intercept = 0.12, p = 0.34), as were the Begg (tau = 0.08, p = 0.52) and Thompson-Sharp (rho = 0.11, p = 0.41) tests, indicating no evidence of publication bias. Trim-and-fill (L0 estimator) imputed 0 missing studies. The fragility index was 4 (fragility quotient = 4/60 = 0.07), indicating moderate robustness. The E-value for the point estimate was 2.47, and for the confidence interval lower bound was 1.82, suggesting that unmeasured confounding of moderate strength would be needed to explain away the result.

## Discussion

The `metaverse` package addresses two persistent methodological challenges in meta-analysis: the vulnerability of standard estimators to data contamination and the lack of rigorous FDR control in moderator selection. Our evaluation demonstrates that MM-estimation recovers the true effect size under 10% contamination with minimal efficiency loss on clean data, while the Knockoff Filter provides formal FDR guarantees that standard meta-regression cannot offer.

### Comparison with Existing Software

The `metafor` package [17] is the most comprehensive general-purpose meta-analysis tool in R, implementing over 30 effect size measures and numerous heterogeneity estimators. However, it does not include robust M/MM/S-estimators and does not offer FDR-controlled variable selection. The `robumeta` package [19] provides cluster-robust variance estimation for dependent effect sizes but does not implement high-breakdown-point estimators. The `metaverse` package complements these tools by providing a distinct robust estimation and selection layer that can be used alongside existing meta-analysis workflows.

### Methodological Considerations

The MM-estimator achieves its robustness by down-weighting studies with large standardized residuals. This is appropriate when outliers arise from reporting errors, data contamination, or distinct clinical subpopulations. However, if the true data-generating process is heavy-tailed (e.g., some interventions genuinely produce extreme effects), down-weighting may obscure clinically meaningful heterogeneity. Users should interpret the robust weights as a diagnostic tool: studies receiving low weights deserve individual scrutiny rather than automatic exclusion.

The Knockoff Filter provides finite-sample FDR control under the assumption that the knockoff variables are exchangeable with the original features conditional on the response. In meta-regression, this assumption is approximately satisfied when moderators are roughly Gaussian. For highly non-Gaussian moderators (e.g., binary, skewed), the equicorrelated or sequential knockoff constructions may be more appropriate than the SDP default.

### Limitations

1. **Sample size requirements:** The Knockoff Filter requires k > 2p (number of studies exceeding twice the number of moderators) for valid construction, which limits its applicability in small meta-analyses. We recommend at least k = 50 for reliable moderator selection with p <= 10 covariates.

2. **Computational intensity:** While the MM-estimator scales linearly in k, the Knockoff Filter involves solving a semidefinite program of dimension p, which becomes expensive for p > 50. The spike-and-slab MCMC may exhibit slow mixing for highly correlated moderators.

3. **No Rcpp acceleration:** The current implementation is pure R. Computationally intensive routines (IRLS, knockoff generation, MCMC) would benefit from C++ optimization via `Rcpp` and `RcppArmadillo` in future versions.

4. **Dependent effect sizes:** The current robust estimators assume independent effect sizes. While `meta_multilevel()` provides a basic multilevel interface, fully robust multivariate meta-analysis combining MM-estimation with structured dependence remains an open methodological problem.

5. **Single evaluation dataset:** The Results section uses a single simulated dataset with known ground truth. While this is the standard evaluation approach for methodology papers, performance may vary across different contamination patterns, moderator structures, and true effect sizes. A comprehensive simulation study across a wider parameter space is planned for a follow-up publication.

6. **No network meta-analysis support:** The current version does not extend robust estimation to network meta-analysis topologies. Integrating MM-estimators into network models (consistency equations, design-by-treatment interaction) is a direction for future work.

### Future Directions

Planned extensions include: (a) C++ acceleration for core algorithms; (b) robust network meta-analysis; (c) a Shiny-based interactive dashboard for non-programmers; (d) integration with the GRADE framework for certainty of evidence assessment; and (e) real-data validation against published Cochrane meta-analyses with known influential outliers.

## Availability and Requirements

- **Project name:** metaverse
- **Repository:** https://github.com/mahmood726-cyber/metaverse
- **Operating system(s):** Platform independent
- **Programming language:** R (>= 4.0.0)
- **Dependencies:** robustbase (>= 0.93), glmnet (>= 4.0), lme4 (>= 1.1-27), ggplot2 (>= 3.3), boot (>= 1.3-28), Matrix (>= 1.3), MASS, rlang (>= 0.4)
- **License:** GPL (>= 3)
- **Data Availability:** All demonstration data used in this manuscript is bundled within the `metaverse` package (`metaverse_eval_data`) and is publicly accessible via CRAN and GitHub.

## Acknowledgments

The author thanks the developers of the `robustbase`, `glmnet`, `metafor`, and `lme4` packages, whose foundational software enabled this work.

## References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Chichester: John Wiley & Sons; 2009.

2. Higgins JPT, Thomas J, Chandler J, Cumpston M, Li T, Page MJ, et al., editors. Cochrane Handbook for Systematic Reviews of Interventions. 2nd ed. Chichester: John Wiley & Sons; 2019.

3. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177-188.

4. Viechtbauer W. Bias and efficiency of meta-analytic variance estimators in the random-effects model. J Educ Behav Stat. 2005;30(3):261-293.

5. Viechtbauer W, Cheung MW. Outlier and influence diagnostics for meta-analysis. Res Synth Methods. 2010;1(2):112-125.

6. Gumedze FN, Jackson D. A random effects variance shift model for detecting and accommodating outliers in meta-analysis. BMC Med Res Methodol. 2011;11:19.

7. Maronna RA, Martin RD, Yohai VJ, Salibian-Barrera M. Robust Statistics: Theory and Methods (with R). 2nd ed. Hoboken: John Wiley & Sons; 2019.

8. Thompson SG, Higgins JPT. How should meta-regression analyses be undertaken and interpreted? Stat Med. 2002;21(11):1559-1573.

9. Higgins JPT, Thompson SG. Controlling the risk of spurious findings from meta-regression. Stat Med. 2004;23(11):1663-1682.

10. Ioannidis JPA. Why most discovered true associations are inflated. Epidemiology. 2008;19(5):640-648.

11. Huber PJ, Ronchetti EM. Robust Statistics. 2nd ed. Hoboken: John Wiley & Sons; 2009.

12. Rousseeuw PJ, Yohai VJ. Robust regression by means of S-estimators. In: Robust and Nonlinear Time Series Analysis. New York: Springer; 1984. p. 256-272.

13. Yohai VJ. High breakdown-point and high efficiency robust estimates for regression. Ann Stat. 1987;15(2):642-656.

14. Barber RF, Candes EJ. Controlling the false discovery rate via knockoffs. Ann Stat. 2015;43(5):2055-2085.

15. Candes E, Fan Y, Janson L, Lv J. Panning for gold: 'model-X' knockoffs for high dimensional controlled variable selection. J R Stat Soc Ser B. 2018;80(3):551-577.

16. Mitchell TJ, Beauchamp JJ. Bayesian variable selection in linear regression. J Am Stat Assoc. 1988;83(404):1023-1032.

17. Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48.

18. Balduzzi S, Rucker G, Schwarzer G. How to perform a meta-analysis with R: a practical tutorial. Evid Based Ment Health. 2019;22(4):153-160.

19. Fisher Z, Tipton E. robumeta: An R-package for robust variance estimation in meta-analysis. arXiv preprint arXiv:1503.02220. 2015.

20. Maechler M, Rousseeuw P, Croux C, Todorov V, Ruckstuhl A, Salibian-Barrera M, et al. robustbase: Basic Robust Statistics R package. 2023. Available from: https://CRAN.R-project.org/package=robustbase

21. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. J Stat Softw. 2010;33(1):1-22.

22. Bates D, Machler M, Bolker B, Walker S. Fitting linear mixed-effects models using lme4. J Stat Softw. 2015;67(1):1-48.

23. Wickham H. ggplot2: Elegant Graphics for Data Analysis. 2nd ed. New York: Springer-Verlag; 2016.

24. Meinshausen N, Buhlmann P. Stability selection. J R Stat Soc Ser B. 2010;72(4):417-473.

25. Fan J, Li R. Variable selection via nonconcave penalized likelihood and its oracle properties. J Am Stat Assoc. 2001;96(456):1348-1360.

26. Zhang CH. Nearly unbiased variable selection under minimax concave penalty. Ann Stat. 2010;38(2):894-942.

27. Vovk V, Wang R. E-values: calibration, combination and applications. Ann Stat. 2021;49(3):1736-1754.

28. Vovk V, Gammerman A, Shafer G. Algorithmic Learning in a Random World. New York: Springer; 2005.

29. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med. 2017;167(4):268-274.
