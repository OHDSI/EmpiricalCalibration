EmpiricalCalibration 3.1.1
==========================

Changes.

1. Making sure we pass R check even if suggested packages are unavailable, as required by CRAN.


EmpiricalCalibration 3.1.0
==========================

Changes

1. Adding the `compareEase()` function to compare EASE of correlated sets of estimates.

2. Adding the `computeCvPoissonRegression()` function to compute critical values for Poisson regressions. 

3. The `plotCalibrationEffect()` function warns if there are values outside the plot limits.

4. Changing the MaxSPRT calibration vignette to use CV calibration instead of LLR calibration.

Bugfixes

1. Fixing `computeTraditionalP()` when using vectors as input.

2. Fix error when no critical value meets requirements (throwing warning instead).


EmpiricalCalibration 3.0.0
==========================

Changes

1. Adding option to make p-values one-sided.

2. Adding empirical calibration of the log likelihood ratio.

3. Adding functions for computing critical values for the log likelihood ratio when performing sequential testing (MaxSPRT).

4. Adding ability to fit null distribution using non-normal approximations of the per-negative control likelihood functions.

5. Faster null distribution fitting using MCMC. Increasing default MCMC iterations for greater stability.

6. Using median instead of mean of posterior distributions when converting MCMC null to error model for greater consistency between calibrated CI and P.


EmpiricalCalibration 2.1.0
==========================

Changes

1. Adding computation of expected absolute systematic error. Can be shown in plots.


EmpiricalCalibration 2.0.2
==========================

Changes

1. Center title on all plots.

Bugfixes

1. Fix build fail in R 4.0.0.


EmpiricalCalibration 2.0.1
==========================

Changes

1. computeTraditionalCi now outputs data frame instead of vector.

Bugfixes

1. convertNullToErrorModel function now adheres to new systematic error model (SD no longer on log scale).


EmpiricalCalibration 2.0.0
==========================

Changes

1. Adding new method for fitting systematic error models that drops the assumption that SD is linear in the log scale, but rather just linear (in the logRr space). Setting that to default, but allowing users to use legacy model if they choose.


EmpiricalCalibration 1.4.1
==========================

Changes

1. Added convertNullToErrorModel function to allow confidence interval calibration using only negative controls (requiring the user to make some assumptions).


EmpiricalCalibration 1.4.0
==========================

Changes

1. Added plot showing effect of confidence interval calibration, similar to p-value plot.

Bugfixes

2. Fixed 'unknown aesthetics' warning when calling plotTrueAndObserved.


EmpiricalCalibration 1.3.6
==========================

Changes

1. Added plots showing expected type 1 error given an estimated empirical null distribution.

2. Closed form solution for RR vs SE plot for faster computation (especially when plotting credible intervals).

Bugfixes

1. Several improvements of the robustness of fitting systematic error models.


EmpiricalCalibration 1.3.1
===========================

Changes

1. Confidence interval calibration model StdDev transformed to log scale to prevent negative StdDev.

2. Confidence interval calibration cross-validation now allows specification of leave-out groups.

3. various new plots for evaluating confidence interval calibration.

4. Added confidence interval calibration vignette.

5. Added example data for confidence interval calibration

Bugfixes

1. None


EmpiricalCalibration 1.2.0
==========================

Changes

1. Ability to add credible intervals to calibration effect plot

2. Plot CI calibration (using leave-one-out cross-validation)

Bugfixes

1. Fixed vignette name in index

2. Removed coverage plot (moved to MethodEvaluation package)


EmpiricalCalibration 1.1.0
==========================

Initial submission to CRAN
