library(testthat)
library(EmpiricalCalibration)

test_that("calibrateP.null argument requirements", {
  
  logRr <- c(0.1,0.3,0.2)
  seLogRr <- c(0.05,0.1)
  null <- fitNull(c(1,2), c(0,0))
  expect_error(calibrateP(null,logRr, seLogRr),
                 regexp = ".*arguments must be of equal length.*")
  })

test_that("calibrateOneP argument requirements", {
  
  logRr <- c(0.1,0.3,0.2)
  seLogRr <- c(0.05,0.1,NA)
  null <- fitNull(c(0,0), c(0,0))
  expect_equal(calibrateP(null,logRr, seLogRr), c(0.045500266,0.002699796,NA))
  expect_equal(calibrateP(null,logRr, seLogRr, twoSided = FALSE), c(0.022750133,0.001349898,NA))
  expect_equal(round(calibrateP(null,logRr, seLogRr, twoSided = FALSE, upper = FALSE),2), c(0.98,1.00,NA))
})


test_that("computeTraditionalP argument options", {
  
  data(sccs)
  positive <- sccs[sccs$groundTruth == 1, ]
  expect_equal(computeTraditionalP(positive$logRr, positive$seLogRr),0)
  expect_equal(computeTraditionalP(positive$logRr, positive$seLogRr, twoSided= FALSE),0)
  expect_equal(computeTraditionalP(positive$logRr, positive$seLogRr, twoSided= FALSE, upper = FALSE),1)
  
})


test_that("fitNullNonNormalLl", {
  
  data(sccs)
  negatives <- sccs[sccs$groundTruth == 0, ]
  null <- fitNull(negatives$logRr, negatives$seLogRr)
  expect_equal(fitNullNonNormalLl(negatives), null)
})

test_that("fitNullNonNormalLl gamma", {
  
  data(sccs)
  negatives <- sccs[sccs$groundTruth == 0, ]
  colnames(negatives)<- c("drugName","mu","gamma", "sigma")
  negatives$mu[1]<- NA
  null <- fitNullNonNormalLl(negatives)
  expect_equivalent(round(null["mean"],1),0)
  expect_equivalent(round(null["sd"],1),0)
  
  colnames(negatives)<- c("drugName","mu","alpha", "sigma")
  null <- fitNullNonNormalLl(negatives)
  expect_equivalent(round(null["mean"],2),0.06)
  expect_equivalent(round(null["sd"],2),0)
  
  colnames(negatives)<- c("drugName","mu","grid", "sigma")
  expect_error(fitNullNonNormalLl(negatives),
               regexp = ".*but not all column names are numeric*" )
  
  colnames(negatives)<- c(1,2,3,4)
  null <- fitNullNonNormalLl(negatives)
  expect_equivalent(round(null["mean"],2),0)
  expect_equivalent(round(null["sd"],2),0.1)
})


