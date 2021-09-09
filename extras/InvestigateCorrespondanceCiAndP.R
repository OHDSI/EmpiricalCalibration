library(EmpiricalCalibration)

# Agreement between P and CI when using MCMC ----------------
simulate <- function(dummy) {
  ncs <- simulateControls() 
  null <- fitMcmcNull(ncs$logRr, ncs$seLogRr, iter = 10000)
  pcs <- simulateControls(trueLogRr = 0.1) 
  
  p <- calibrateP(null, pcs$logRr, pcs$seLogRr)
  model <- convertNullToErrorModel(null)
  # chain <- attr(null, "mcmc")$chain
  # # model <- c(mean(chain[, 1]), 1, 1/sqrt(mean(chain[, 2])), 0) #36 errors
  # model <- c(median(chain[, 1]), 1, 1/sqrt(median(chain[, 2])), 0) # 3 errors
  # names(model) <- c("meanIntercept", "meanSlope", "sdIntercept", "sdSlope")
  # class(model) <- "systematicErrorModel"
  
  ci <- calibrateConfidenceInterval(pcs$logRr, pcs$seLogRr, model)
  
  significantP <- p$p < 0.05
  significantCi <- ci$logLb95Rr > 0 | ci$logUb95Rr < 0
  errors <- sum(significantP & !significantCi) + sum(!significantP & significantCi)
  return(errors)
}
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "EmpiricalCalibration")
errors <- ParallelLogger::clusterApply(cluster, 1:100, simulate)
errors <- do.call(c, errors)
sum(errors)

# Agreement between P using MCMC and P using asymptotics ----------------
simulate <- function(dummy) {
  ncs <- simulateControls() 
  null <- fitMcmcNull(ncs$logRr, ncs$seLogRr, iter = 10000)
  pcs <- simulateControls(trueLogRr = 0.1) 
  
  p <- calibrateP(null, pcs$logRr, pcs$seLogRr)
  
  null2 <- c(null[1], 1/sqrt(null[2]))
  class(null2) <- "null"
  
  # chain <- attr(null, "mcmc")$chain
  # null2 <- c(mean(chain[, 1]), 1/sqrt(mean(chain[, 2]))) #33 errors
  # # null2 <- c(median(chain[, 1]), 1/sqrt(median(chain[, 2]))) # 1 error
  # class(null2) <- "null"
  
  p2 <- calibrateP(null2, pcs$logRr, pcs$seLogRr)
  
  significantP <- p$p < 0.05
  significantP2 <- p2 < 0.05
  errors <- sum(significantP & !significantP2) + sum(!significantP & significantP2)
  return(errors)
}
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "EmpiricalCalibration")
errors <- ParallelLogger::clusterApply(cluster, 1:100, simulate)
errors <- do.call(c, errors)
sum(errors)
