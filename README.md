EmpiricalCalibration
====================

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/EmpiricalCalibration)](http://cran.r-project.org/package=EmpiricalCalibration)

Introduction
============

This R package contains routines for performing empirical calibration of observational study estimates. By using a set of negative control hypotheses we can estimate the empirical null distribution of a particular observational study setup. This empirical null distribution can be used to compute a calibrated p-value, which reflects the probability of observing an estimated effect size when the null hypothesis is true taking both random and systematic error into account, as described in the paper [Interpreting observational studies: why empirical calibration is needed to correct p-values.] (http://dx.doi.org/10.1002/sim.5925).

Features
========
- Estimate the empirical null distribution given the effect estimates of a set of negative controls 
- Estimate the calibrated p-value of a given hypothesis given the  estimated empirical null distribution
- Produce various plots for evaluating the empirical calibration
- Contains the data sets from the paper for illustration

Screenshots and examples
========================
<img src="https://github.com/OHDSI/EmpiricalCalibration/blob/master/extras/plot.png" alt="Calibration effect plot" title="Calibration effect plot" />

```r
data(sccs) #Load one of the included data sets
negatives <- sccs[sccs$groundTruth == 0,] #Select the negative controls
null <- fitNull(negatives$logRr,negatives$seLogRr) #Fit the null distribution
positive <- sccs[sccs$groundTruth == 1,]  #Select the positive control

#Create the plot above:
plotCalibrationEffect(negatives$logRr,negatives$seLogRr,positive$logRr,positive$seLogRr,null)

#Compute the calibrated p-value:
calibrateP(positive$logRr,positive$seLogRr, null) #Compute calibrated p-value
[1] 0.8390598
```

Technology
==========
This is a pure R package.

System requirements
===================
Requires [R](http://cran.r-project.org/) (version 3.1.0 or newer).

Getting Started
===============
In R, use the following commands to install the latest stable version from CRAN:

```r
install.packages("EmpiricalCalibration")
```

To install the latest development version directly from GitHub, use:

```r
install.packages("devtools")
library(devtools)
install_github("ohdsi/EmpiricalCalibration")
```
  
Getting Involved
================
* Vignette: [Empirical calibration of p-values](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/master/inst/doc/EmpiricalPCalibrationVignette.pdf)
* Vignette: [Empirical calibration of confidence intervals](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/master/inst/doc/EmpiricalCiCalibrationVignette.pdf)
* Package manual: [EmpiricalCalibration.pdf](https://raw.githubusercontent.com/OHDSI/EmpiricalCalibration/master/extras/EmpiricalCalibration.pdf) 
* Developer questions/comments/feedback: <a href="http://forums.ohdsi.org/c/developers">OHDSI Forum</a>
* We use the <a href="../../issues">GitHub issue tracker</a> for all bugs/issues/enhancements
  
License
=======
EmpiricalCalibration is licensed under Apache License 2.0

Development
===========
This package has been developed in RStudio.

### Development status
[![Build Status](https://travis-ci.org/OHDSI/EmpiricalCalibration.svg?branch=master)](https://travis-ci.org/OHDSI/EmpiricalCalibration)

This package is ready for use.

Acknowledgements
================
Martijn Schuemie is the author of this package.
