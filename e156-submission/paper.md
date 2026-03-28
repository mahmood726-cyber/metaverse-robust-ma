Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

metaverse: Robust Meta-Analysis with Outlier-Resistant Estimation and Variable Selection

Can robust estimation methods protect meta-analytic conclusions from outlier contamination while maintaining valid inference under standard conditions? We developed the metaverse R package implementing M, MM, S, and tau-scale estimators with contamination models alongside knockoff filters, spike-and-slab selection, penalized regression, conformal prediction, and post-selection inference for moderators. The package provides 12 plot types including forest, funnel, radial, and influence diagnostics, plus effect-size converters covering 10 metrics and power analysis for study planning. Under 10 percent contamination across 60 studies the pooled SMD bias was 0.02, 95% CI 0.00-0.04, for the MM-estimator versus 0.15 for DerSimonian-Laird. Knockoff-filtered variable selection maintained FDR below the nominal 10 percent threshold while correctly detecting both true moderators in the evaluation dataset. The framework integrates with metafor workflows and supports dependent effect sizes through multilevel modeling via the lme4 package. A limitation is that robust estimators require approximately 20 studies to reliably outperform standard methods under low contamination.

Outside Notes

Type: methods
Primary estimand: Pooled SMD bias
App: metaverse R package v0.1.0
Data: Simulated 60-study dataset with 10% contamination
Code: https://github.com/mahmood726-cyber/metaverse-robust-ma
Version: 0.1.0
Validation: DRAFT

References

1. Van den Noortgate W, Lopez-Lopez JA, Marin-Martinez F, Sanchez-Meca J. Three-level meta-analysis of dependent effect sizes. Behav Res Methods. 2013;45:576-594.
2. Assink M, Wibbelink CJM. Fitting three-level meta-analytic models in R: a step-by-step tutorial. Quant Methods Psychol. 2016;12(3):154-174.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
