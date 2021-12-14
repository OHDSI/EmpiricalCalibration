# Some code to explore what we can do if we have estimates from different
# analyses (e.g cohort method and SCCS) for multiple time-at-risk (TAR) windows.

library(Cyclops)
library(ParallelLogger)

settings <- list(
  nTarget = 1000,              # Number of subjects in the target population
  nComparator = 10000,         # Number of subjects in the comparator (counterfactual) population
  backgroundRate = 0.001,      # Poisson background rate of the outcome    
  rr = 2,                      # Relative risk during the true TAR when exposed.
  trueTarStart = 0,            # Start of the true TAR (relative to exposure)
  trueTarEnd = 21,             # End of the true TAR (relative to exposure)
  tarStart = c(0, 0, 0, 0),    # Start times for the various TARs used in the analyses
  tarEnd = c(7, 21, 42, 90),   # End times for the various TARs used in the analyses
  nMethods = 2                 # Number of different methods (here simulated as different sampling of the comparator)
)

# Simulation ----------------------------------------------------------------

runSimulation <- function(iteration, settings) {
  
  computeEventsPerTar <- function(n, rr) {
    cutOffPoints <- c(settings$trueTarStart, settings$trueTarEnd, settings$tarStart, settings$tarEnd)
    cutOffPoints <- sort(unique(cutOffPoints))
    intervals <- data.frame(start = head(cutOffPoints, -1),
                            end = tail(cutOffPoints, -1))
    intervals$trueTar <- intervals$start >= settings$trueTarStart & 
      intervals$end <= settings$trueTarEnd 
    
    intervals$lambda <- n * ifelse(intervals$trueTar, rr * settings$backgroundRate, settings$backgroundRate) * (intervals$end - intervals$start + 1)
    intervals$events <- rpois(nrow(intervals), intervals$lambda)
    
    tarEvents <- sapply(1:length(settings$tarStart), function(i) sum(intervals$events[intervals$start >= settings$tarStart[i] & intervals$end <= settings$tarEnd[i]]))
    return(tarEvents)
  }
  
  computeEstimates <- function(method, targetEvents) {
    timePerPerson <- settings$tarEnd - settings$tarStart + 1
    estimates <- data.frame(tar = 1:length(settings$tarStart),
                            start = settings$tarStart,
                            end =settings$ tarEnd,
                            targetEvents = targetEvents,
                            comparatorEvents = computeEventsPerTar(settings$nComparator, 1),
                            targetTime = timePerPerson * settings$nTarget,
                            comparatorTime = timePerPerson * settings$nComparator,
                            logRr = 0,
                            logCi95Lb = 0,
                            logCi95Ub = 0,
                            seLogRr = 0)
    for (i in 1:length(settings$tarStart)) {
      x <- c(1, 0)
      y <- c(estimates$targetEvents[i], estimates$comparatorEvents[i])
      logTime <- log(c(estimates$targetTime[i], estimates$comparatorTime[i]))
      
      cyclopsData <- createCyclopsData(y ~ x + offset(logTime), modelType = "pr")
      fit <- fitCyclopsModel(cyclopsData)
      
      # Simplistic approach: assume likelihood is normally distributed, so can be 
      # expressed as logRr and seLogRr:
      logRr <- coef(fit)["x"]
      logCi <- confint(fit, parm = "x")
      estimates$logRr[i] <- logRr
      estimates$logCi95Lb[i] <- logCi[2]
      estimates$logCi95Ub[i] <- logCi[3]
      estimates$seLogRr[i] <- (logCi[3] - logCi[2])/(2 * qnorm(0.975))
    }
    estimates$method <- method
    return(estimates)
  }
  
  targetEvents <- computeEventsPerTar(settings$nTarget, settings$rr)
  estimates <- lapply(1:settings$nMethods, computeEstimates, targetEvents = targetEvents)   
  estimates <- do.call(rbind, estimates) 
  
  computeAverage <- function(models) {
    # Simplistic approach to averaging models: non-Bayesian averaging, assuming 
    # the average is still convex (which will be untrue in many scenarios):
    
    ll <- function(x, models) {
      ll <- -sapply(x, function(x) sum(dnorm(x, mean = models$logRr, sd = models$seLogRr, log = TRUE)) / settings$nMethods)
      return(ll)
    } 
    a <- 0.05
    fit <- nlm(ll, mean(models$logRr), models = models)
    logRr <- fit$estimate
    threshold <- fit$minimum - qchisq(1 - a, df = 1)/2
    
    precision <- 1e-07
    
    # Binary search for upper bound
    L <- logRr
    H <- 10
    ub <- Inf
    while (H >= L) {
      M <- L + (H - L)/2
      llM <- -ll(M, models)
      metric <- threshold - llM
      if (metric > precision) {
        H <- M
      } else if (-metric > precision) {
        L <- M
      } else {
        ub <- M
        break
      }
      if (M == logRr) {
        warn("Error finding upper bound")
        break
      } else if (M == 10) {
        warn("Confidence interval upper bound out of range")
        break
      }
    }
    
    # Binary search for lower bound
    L <- -10
    H <- logRr
    lb <- -Inf
    while (H >= L) {
      M <- L + (H - L)/2
      llM <- -ll(M, models)
      metric <- threshold - llM
      if (metric > precision) {
        L <- M
      } else if (-metric > precision) {
        H <- M
      } else {
        lb <- M
        break
      }
      if (M == logRr) {
        warn("Error finding lower bound")
        break
      } else if (M == -10) {
        warn("Confidence interval lower bound out of range")
        break
      }
    }
    result <- data.frame(logRr = logRr,
                         logCi95Lb = lb,
                         logCi95Ub = ub,
                         seLogRr = (ub - lb)/(2 * qnorm(0.975)))
    return(result)
  }
  
  averageModels <- data.frame()
  for (i in 1:length(settings$tarStart)) {
    models <- estimates[estimates$tar == i, ]
    averageModel <- computeAverage(models)
    averageModel$tar <- i
    averageModel$start <- models$start[1]
    averageModel$end <- models$end[1]
    averageModel$method <- "average"
    averageModels <- rbind(averageModels, averageModel)
  }
  result <- rbind(estimates[, c("tar", "start", "end", "logRr", "logCi95Lb", "logCi95Ub", "seLogRr", "method")],
                  averageModels[, c("tar", "start", "end", "logRr", "logCi95Lb", "logCi95Ub", "seLogRr", "method")])
  result$iteration <- iteration
  return(result)
}

cluster <- makeCluster(10)
clusterRequire(cluster, "Cyclops")
results <- clusterApply(cluster, 1:1000, runSimulation, settings = settings)
stopCluster(cluster)
results <- do.call(rbind, results)

# Compute coverage of the 95% CI per method (inc. model averaging):
for (tar in 1:length(settings$tarStart)) {
  for (method in unique(results$method)) {
    subset <- results[results$method == method & results$tar == tar, ]
    coverage <- mean(subset$logCi95Lb <= log(settings$rr) & subset$logCi95Ub >= log(settings$rr))
    writeLines(sprintf("TAR: %s-%s, method: %s, coverage: %0.3f", 
                       settings$tarStart[tar],
                       settings$tarEnd[tar],
                       method,
                       coverage))
  }
}

