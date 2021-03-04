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
