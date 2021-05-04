library(survival)

# Cox proportional hazards regression (cohort method) -----------------------------------------------------
parameters <- data.frame(n = 100000, # Number of subjects
                         pExposure = 0.1, # Probability of being in target cohort
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
      return(data.frame(logRr = NA,
                        lb = 0,
                        llr = 0,
                        events = sum(data$outcome),
                        exposedTime = sum(data$time[data$exposure]),
                        unexposedTime = sum(data$time[!data$exposure])))
    } else {
      llAproximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit = fit,
                                                                 parameter = "exposureTRUE",
                                                                 approximation = "custom")
                                               
      
      if (useCalibration) {
        null <- c(parameters$nullMu, parameters$nullSigma)
        names(null) <- c("mean", "sd")
        class(null) <- "null"
        llr <- EmpiricalCalibration::calibrateLlr(null, llAproximation, twoSided = FALSE, upper = TRUE)
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
      return(data.frame(llr,
                        events = sum(data$outcome),
                        exposedTime = sum(data$time[data$exposure]),
                        unexposedTime = sum(data$time[!data$exposure])))
    }
  }     
  tIndex <- runif(parameters$n,  0, parameters$maxT)
  exposure <- runif(parameters$n) < parameters$pExposure
  systematicError <- rnorm(n = 1, mean = parameters$nullMu, sd = parameters$nullSigma)
  tOutcome <- rexp(parameters$n,  parameters$backgroundHazard * (1 + ((parameters$rr - 1) * exp(systematicError) * exposure)))
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
cluster <- ParallelLogger::makeCluster(2)
ParallelLogger::clusterRequire(cluster, "survival")

mean(unlist(ParallelLogger::clusterApply(cluster, 1:1000, simulate, parameters = parameters, useCalibration = TRUE, looks = 10)), na.rm = TRUE)

ParallelLogger::stopCluster(cluster)


