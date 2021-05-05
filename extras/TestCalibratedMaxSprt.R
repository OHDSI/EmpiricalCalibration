library(survival)

# Cox proportional hazards regression (cohort method) -----------------------------------------------------
parameters <- data.frame(n = 100000, # Number of subjects
                         pExposure = 0.5, # Probability of being in target cohort
                         backgroundHazard = 0.0001,
                         tar = 10, # Time at risk for each exposure
                         rr = 1, # Relative risk (hazard ratio)
                         nullMu = 0.2, # Null distribution mean (at log HR scale)
                         nullSigma = 0.2, # Null distribution SD (at log HR scale)
                         maxT = 100) 

simulate <- function(seed, parameters, useCalibration = TRUE, looks = 10) {
  set.seed(seed)
  
  computeAtT <- function(t) {
    truncatedTime <- time
    idxTruncated <- tIndex + time > t
    truncatedTime[idxTruncated] <- t - tIndex[idxTruncated]
    truncatedOutcome <- outcome
    truncatedOutcome[idxTruncated] <- 0
    data <- data.frame(time = truncatedTime,
                       outcome = truncatedOutcome,
                       exposure = exposure)
    data <- data[data$time > 0, ]
    cyclopsData <- Cyclops::createCyclopsData(Surv(time, outcome) ~ exposure , modelType = "cox", data = data)
    fit <- Cyclops::fitCyclopsModel(cyclopsData, control = Cyclops::createControl(seed = seed))
    if (fit$return_flag != "SUCCESS") {
      return(data.frame(llr = 0,
                        events = sum(data$outcome),
                        exposedTime = sum(data$time[data$exposure]),
                        unexposedTime = sum(data$time[!data$exposure])))
    } else {
      if (useCalibration) {
        llApproximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit = fit,
                                                                   parameter = "exposureTRUE",
                                                                   approximation = "grid")
        null <- c(parameters$nullMu, parameters$nullSigma)
        names(null) <- c("mean", "sd")
        class(null) <- "null"
      
        llr <- EmpiricalCalibration::calibrateLlr(null = null, 
                                                  likelihoodApproximation = llApproximation, 
                                                  twoSided = FALSE, 
                                                  upper = TRUE)
      } else {
        if (coef(fit)["exposureTRUE"] > 0) {
          llNull <- Cyclops::getCyclopsProfileLogLikelihood(object = fit,
                                                            parm = "exposureTRUE",
                                                            x = 0,
                                                            includePenalty = FALSE)$value
          llr <- fit$log_likelihood - llNull
        } else {
          llr <- 0
        }
      }
    }
    return(data.frame(llr,
                      events = sum(data$outcome),
                      exposedTime = sum(data$time[data$exposure]),
                      unexposedTime = sum(data$time[!data$exposure])))
  }     
  tIndex <- runif(parameters$n,  0, parameters$maxT)
  exposure <- runif(parameters$n) < parameters$pExposure
  systematicError <- rnorm(n = 1, mean = parameters$nullMu, sd = parameters$nullSigma)
  hazard <- ifelse(exposure, parameters$backgroundHazard * parameters$rr * exp(systematicError), parameters$backgroundHazard)
  tOutcome <- rexp(parameters$n,  hazard)
  outcome <- tOutcome < parameters$tar
  time <- rep(parameters$tar, parameters$n)
  time[outcome] <- tOutcome[outcome]
  t <- seq(0, parameters$maxT, length.out = looks + 1)[-1]                 
  results <- purrr::map_dfr(t, computeAtT)
  
  # Perform sequential testing
  sampleSizeUpperLimit <- max(results$events, na.rm = TRUE)
  if (sampleSizeUpperLimit <= 5) {
    return(FALSE)
  }
  events <- results$events
  if (looks > 1) {
    events[2:looks] <- events[2:looks] - events[1:(looks - 1)]
    events <- events[events != 0]
  }
  cv <- Sequential::CV.Binomial(N = sampleSizeUpperLimit,
                                M = 1,
                                z = max(results$unexposedTime) / max(results$exposedTime),
                                GroupSizes = events)$cv
  return(any(results$llr > cv, na.rm = TRUE))
}

# Compute type I error (probability of a signal when the null is true). Should be 0.05:
cluster <- ParallelLogger::makeCluster(20)
ParallelLogger::clusterRequire(cluster, "survival")


# Scenario 1: no systematic error, 1 look
parameters$nullMu <- 0
parameters$nullSigma <- 0
mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 1)), na.rm = TRUE)
# [1] 0.05

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = FALSE, looks = 1)), na.rm = TRUE)
# [1] 0.055


# Scenario 2: no systematic error, 10 sequential looks
mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 10)), na.rm = TRUE)
# [1] 0.048

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = FALSE, looks = 10)), na.rm = TRUE)
# [1] 0.064


# Scenario 3: systematic error, 1 look
parameters$nullMu <- 0.2
parameters$nullSigma <- 0.2
mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 1)), na.rm = TRUE)
# [1] 0.045

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = FALSE, looks = 1)), na.rm = TRUE)
# [1] 0.324


# Scenario 4: systematic error, 10 sequential looks
mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 10)), na.rm = TRUE)
# [1] 0.04

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = FALSE, looks = 10)), na.rm = TRUE)
# [1] 0.258


# Scenario 5: no systematic error, 1 look, imbalanced exposure cohorts
parameters$nullMu <- 0
parameters$nullSigma <- 0
parameters$pExposure <- 0.1
mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 1)), na.rm = TRUE)
# [1] 0.04

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = FALSE, looks = 1)), na.rm = TRUE)
# [1] 0.041

# Scenario 6: systematic error, 10 sequential looks, imbalanced exposure cohorts
parameters$nullMu <- 0.2
parameters$nullSigma <- 0.2
parameters$pExposure <- 0.1
mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 10)), na.rm = TRUE)
# [1] 0.043

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = FALSE, looks = 10)), na.rm = TRUE)
# [1] 0.146

ParallelLogger::stopCluster(cluster)



# Abstract Poisson (historic comparator method) -------------------------------------------------------------
parameters <- data.frame(n = 1000000, # Subjects
                         rate = 0.0001,
                         nullMu = 0.2, # Null distribution mean (at log IRR scale)
                         nullSigma = 0.2) # Null distribution SD (at log IRR scale))

simulate <- function(seed, parameters, useCalibration = TRUE, cv) {
  set.seed(seed)
  
  llr <- function(observed, expected) {
    if (observed <= expected) {
      return(0)
    } else {
      if (useCalibration) {
        bounds = c(log(0.1), log(10))
        x <- seq(bounds[1], bounds[2], length.out = 1000)
        ll <- dpois(observed, expected * exp(x), log = TRUE)
        names(ll) <- x
        null <- c(parameters$nullMu, parameters$nullSigma)
        names(null) <- c("mean", "sd")
        class(null) <- "null"
        return(EmpiricalCalibration::calibrateLlr(null = null, likelihoodApproximation = ll))
      } else {
        return(dpois(observed, observed, log = TRUE) - dpois(observed, expected, log = TRUE))
      }
    }
  }
  systematicError <- rnorm(n = 1, mean = parameters$nullMu, sd = parameters$nullSigma)
  observed <- sum(rpois(parameters$n, parameters$rate * exp(systematicError)))
  expected <- parameters$n * parameters$rate
  llr <- llr(observed = observed, expected = expected)
  return(llr > cv)
}

# Compute type I error (probability of a signal when the null is true). Should be 0.05:
expected <- parameters$n * parameters$rate
cv <- Sequential::CV.Poisson(SampleSize = expected,
                             alpha = 0.05,
                             M = 1,
                             GroupSizes = c(expected))

cluster <- ParallelLogger::makeCluster(10)
mean(unlist(ParallelLogger::clusterApply(cluster, 1:10000, simulate, parameters = parameters, cv = cv, useCalibration = TRUE)), na.rm = TRUE)
# [1] 0.0471

mean(unlist(ParallelLogger::clusterApply(cluster, 1:10000, simulate, parameters = parameters, cv = cv, useCalibration = FALSE)), na.rm = TRUE)
# [1] 0.5726

ParallelLogger::stopCluster(cluster)
