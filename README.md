EmpiricalCalibration
====================

[![Build Status](https://github.com/OHDSI/EmpiricalCalibration/workflows/R-CMD-check/badge.svg)](https://github.com/OHDSI/EmpiricalCalibration/actions?query=workflow%3AR-CMD-check)
[![codecov.io](https://codecov.io/github/OHDSI/EmpiricalCalibration/coverage.svg?branch=main)](https://app.codecov.io/github/OHDSI/EmpiricalCalibration?branch=main)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/EmpiricalCalibration)](https://cran.r-project.org/package=EmpiricalCalibration)
[![CRAN_Status_Badge](http://cranlogs.r-pkg.org/badges/EmpiricalCalibration)](https://cran.r-project.org/package=EmpiricalCalibration)

EmpiricalCalibration is part of [HADES](https://ohdsi.github.io/Hades/).

Introduction
============

This R package contains routines for performing empirical calibration of observational study estimates. By using a set of negative control hypotheses we can estimate the empirical null distribution of a particular observational study setup. This empirical null distribution can be used to compute a calibrated p-value, which reflects the probability of observing an estimated effect size when the null hypothesis is true taking both random and systematic error into account, as described in the paper [Interpreting observational studies: why empirical calibration is needed to correct p-values](http://dx.doi.org/10.1002/sim.5925). 

Also supported is empirical calibration of confidence intervals, based on the results for a set of negative and positive controls, as described in the paper [Empirical confidence interval calibration for population-level effect estimation studies in observational healthcare data](https://doi.org/10.1073/pnas.1708282114).

Features
========
- Estimate the empirical null distribution given the effect estimates of a set of negative controls. 
- Estimate the calibrated p-value of a given hypothesis given the estimated empirical null distribution.
- Estimate a systematic error distribution given the effect estimates for a set of negative and positive controls.
- Estimate the calibrated confidence interval for a given estimate given the systematic error distribution.
- Estimate a calibrated log likelihood ratio, for use in maximum sequential probability ratio testing (MaxSPRT).
- Produce various plots for evaluating the empirical calibration.
- Contains the data sets from the papers for illustration.

Screenshots and examples
========================
<img src="https://github.com/OHDSI/EmpiricalCalibration/raw/main/extras/plot.png" alt="Calibration effect plot" title="Calibration effect plot" />

```r
data(sccs) #Load one of the included data sets
negatives <- sccs[sccs$groundTruth == 0,] #Select the negative controls
null <- fitNull(logRr = negatives$logRr, seLogRr = negatives$seLogRr) #Fit the null distribution
positive <- sccs[sccs$groundTruth == 1,]  #Select the positive control

#Create the plot above:
plotCalibrationEffect(logRrNegatives = negatives$logRr,
                      seLogRrNegatives = negatives$seLogRr,
                      logRrPositives = positive$logRr,
                      seLogRrPositives = positive$seLogRr,
                      null = null)

#Compute the calibrated p-value:
calibrateP(null = null, logRr = positive$logRr, seLogRr = positive$seLogRr) #Compute calibrated p-value
[1] 0.8390598
```

Technology
==========
This is a pure R package.

System requirements
===================
Requires R (version 3.1.0 or newer).

Installation
============
In R, use the following commands to install the latest stable version from CRAN:

```r
install.packages("EmpiricalCalibration")
```

To install the latest development version directly from GitHub, use:

```r
install.packages("remotes")
library(remotes)
install_github("ohdsi/EmpiricalCalibration", ref = "develop")
```
  
User Documentation
==================
Documentation can be found on the [package website](https://ohdsi.github.io/EmpiricalCalibration/).

PDF versions of the documentation is also available:

* Vignette: [Empirical calibration of p-values](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/main/inst/doc/EmpiricalPCalibrationVignette.pdf)
* Vignette: [Empirical calibration of confidence intervals](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/main/inst/doc/EmpiricalCiCalibrationVignette.pdf)
* Vignette: [Empirical calibration and MaxSPRT](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/main/inst/doc/EmpiricalMaxSprtCalibrationVignette.pdf)
* Package manual: [EmpiricalCalibration.pdf](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/main/extras/EmpiricalCalibration.pdf) 

Support
=======
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="https://github.com/OHDSI/EmpiricalCalibration/issues">GitHub issue tracker</a> for all bugs/issues/enhancements

Contributing
============
Read [here](https://ohdsi.github.io/Hades/contribute.html) how you can contribute to this package.

License
=======
EmpiricalCalibration is licensed under Apache License 2.0

Development
===========
This package has been developed in RStudio.

### Development status

This package is ready for use.
