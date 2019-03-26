EmpiricalCalibration 1.4.1
==========================

NEW FEATURES

* Added convertNullToErrorModel function to allow confidence interval calibration using only negative controls (requiring the user to make some assumptions).

EmpiricalCalibration 1.4.0
==========================

NEW FEATURES

* Added plot showing effec of confidence interval calibration, similar to p-value plot.

BUG FIXES

* Fixed 'unknown aesthetics' warning when calling plotTrueAndObserved.

EmpiricalCalibration 1.3.6
==========================

NEW FEATURES

* Added plots showing expected type 1 error given an estimated empirical null distribution.

* Closed form solution for RR vs SE plot for faster computation (especially when plotting credible intervals).

BUG FIXES

* Several improvements of the robustness of fitting systematic error models.

EmpiricalCalibration 1.3.1
===========================

NEW FEATURES

* Confidence interval calibration model StdDev transformed to log scale to prevent negative StdDev.

* Confidence interval calibration cross-validation now allows specification of leave-out groups.

* various new plots for evaluating confidence interval calibration.

* Added confidence interval calibration vignette.

* Added example data for confidence interval calibration

BUG FIXES

* None


EmpiricalCalibration 1.2.0
==========================

NEW FEATURES

* Ability to add credible intervals to calibration effect plot

* Plot CI calibration (using leave-one-out cross-validation)

BUG FIXES

* Fixed vignette name in index

* Removed coverage plot (moved to MethodEvaluation package)


EmpiricalCalibration 1.1.0
==========================

Changes: initial submission to CRAN
