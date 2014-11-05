EmpiricalCalibration
====================

This R package contains routines for performing empirical calibration of observational study estimates. By using a set of negative control hypotheses we can estimate the empirical null distribution of a particular observational study setup. This empirical null distribution can be used to compute a calibrated p-value, which reflects the probability of observing an estimated effect size when the null hypothesis is true taking both random and systematic error into account, as described in the paper [Interpreting observational studies: why empirical calibration is needed to correct p-values.] (http://dx.doi.org/10.1002/sim.5925).

Features
========

- Estimate the empirical null distribution given the effect estimates of a set of negative controls 
- Estimate the calibrated p-value of a given hypothesis given the  estimated empirical null distribution
- Produce various plots for evaluating the empirical calibration
- Contains the data sets from the paper for illustration

Requirements
============
Requires [R](http://cran.r-project.org/) with the packages ggplot2 and MASS installed. If you want to perform empirical calibration of your estimates, you need to have a set of negative control hypothesis for which the effect sizes and standard errors have also been estimated.

Screenshots and examples
========================
<img src="https://github.com/OHDSI/EmpiricalCalibration/blob/master/man/plot.png" alt="Calibration effect plot" title="Calibration effect plot" />

```r
data(sccs) #Load one of the included data sets
negatives <- sccs[sccs$groundTruth == 0,] #Select the negative controls
null <- fitNull(negatives$logRr,negatives$seLogRr) #Fit the null distribution
positive <- sccs[sccs$groundTruth == 1,]  #Select the positive control / hypothesis of interest
calibrateP(positive$logRr,positive$seLogRr, null) #Compute calibrated p-value
[1] 0.8390598
```

Getting Started
===============
In R, use the following commands to install this package:

```r
install.packages("devtools")
library(devtools)
install_github("ohdsi/EmpiricalCalibration")
```
  
Documentation
=============
Once installed, use

```r
library(EmpiricalCalibration)
?EmpiricalCalibration
```
to get an overview of the functions and dataset included in this packages.
  
License
=======
Achilles is licensed under Apache License 2.0

Development
===========
This package has been developed in RStudio.
###Development status
This package is ready for use.
