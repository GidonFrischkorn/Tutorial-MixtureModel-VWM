# Tutorial: Mixture Models for visual working memory tasks

This repository contains R scripts &amp; raw data to illustrate how to estimate widely used mixture models (e.g. the Zhang &amp; Luck (2008) or the Bays et al. (2009) mixture models) for visual working memory tasks.

To ensure that all relative paths work, we have set up an R project with the following file structure:

-> Tutorial-MixtureModel-VWM.Rproj: open this project file in RStudio. That way all relative paths work 
- scripts: folder containing all R scripts
- output: folder for output files of R scripts (e.g. RData files containing brms model fits)
- figures: folder to save figures from R scripts
- data: folder containing benchmark data sets to evaluate mixture models

Please note: Some of the output files exceeded the size limit for files allowed by GitHub, these are avaliable via: https://osf.io/vsrz4/

## Package versions

The model fits on OSF were made with development versions of bmm from 0.5.3.9000 to 1.0.1.9000 and with brms 2.21 to 2.22. On 10 November 2026 we refit the IMMbsc example with bmm 1.0.1, because bmm 1.0 changed the parametrisation of the IMM and the stored fit predated that change. If you downloaded `fit_E5_OL2017_IMMbsc.rds` or `E5_bridge_imm_bsc.rds` before that date, delete both files and download them again (`LoadResultsFiles.R` skips files that already exist).

Newer bmm versions may set different default priors than the ones stored with the fits, for example on group-level SDs. Results then differ slightly; `prior_summary(fit)` on a stored fit shows which priors were used.

With bmm 1.4.1 and the current development version, Example 4 with monotonic constraints fails at the first iteration (popov-lab/bmm#694). bmm 1.0.1 runs it.
