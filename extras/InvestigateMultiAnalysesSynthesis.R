# Some code to explore what we can do if we have estimates from different
# analyses (e.g cohort method and SCCS) for multiple time-at-risk (TAR) windows.

library(Cyclops)
library(ParallelLogger)
library(dplyr)
library(tidyr)

settings <- list(
  nTarget = 1000,              # Number of subjects in the target population
  nComparator = 10000,         # Number of subjects in the comparator (counterfactual) population
  backgroundRate = 0.001,      # Poisson background rate of the outcome    
  rr = 2,                      # Relative risk during the true TAR when exposed
  tar = 21,                    # Time-at-risk per subject
  nMethods = 2,                # Number of different methods (here simulated as different sampling of the comparator)
  nBootstrap = 100             # Number of samples for the bootstrap
)

# Simulation ----------------------------------------------------------------

runSimulation <- function(iteration, settings) {
  
  simulatePersons <- function(groupId, n, rr) {
    persons <- tibble(events = rpois(n, settings$backgroundRate * rr * settings$tar),
                      groupId = groupId)
    return(persons)
  }
  
  computeEstimatePerMethod <- function(method, persons) {
    subset <- persons %>%
      filter(groupId == 0 | groupId == method) %>%
      mutate(x = groupId == 0,
             logTime = log(settings$tar))
    cyclopsData <- createCyclopsData(events ~ x + offset(logTime), data = subset, modelType = "pr")
    fit <- fitCyclopsModel(cyclopsData)
    
    # Simplistic approach: assume likelihood is normally distributed, so can be 
    # expressed as logRr and seLogRr:
    logRr <- coef(fit)["xTRUE"]
    logCi <- confint(fit, parm = "xTRUE")
    estimate <- tibble(
      logRr = logRr,
      logCi95Lb = logCi[2],
      logCi95Ub = logCi[3],
      seLogRr = (logCi[3] - logCi[2])/(2 * qnorm(0.975)),
      method = method
    )
    return(estimate)
  }
  
  computeEstimates <- function(persons) {
    estimates <- lapply(1:settings$nMethods, computeEstimatePerMethod, persons = persons)
    estimates <- bind_rows(estimates)
    return(estimates)
  }
  
  performBootstrap <- function(sampleId, persons) {
    sampledPersons <- persons[sample.int(n = nrow(persons), replace = TRUE), ]
    estimates <- computeEstimates(sampledPersons) %>%
      mutate(sampleId = sampleId)
    return(estimates)
  }
  
  computeCorrelations <- function(boostrapEstimates) {
    # For efficiency: pivot so we have one column per method:
    pivotedTable <- boostrapEstimates %>%
      pivot_wider(names_from = method, names_prefix = "method_", id_cols = "sampleId", values_from = "logRr")
    # Compute correlation for all pairs of methods:
    correlations <- expand.grid(methodA = 1:settings$nMethods, methodB = 1:settings$nMethods) %>%
      as_tibble() %>%
      filter(methodA < methodB) 
    correlations$correlation <- sapply(1:nrow(correlations), function(i) cor(pivotedTable[, sprintf("method_%d", correlations$methodA[i])], 
                                                                             pivotedTable[, sprintf("method_%d", correlations$methodB[i])]))
    return(correlations)
  }
  
  # Creating one big population with target (groupId = 0) and all comparators:
  targetPersons <- simulatePersons(groupId = 0, n = settings$nTarget, rr = settings$rr)
  comparatorPersons <- lapply(1:settings$nMethods, simulatePersons, n = settings$nComparator, rr = 1)
  persons <- bind_rows(c(list(targetPersons), comparatorPersons))
  
  # Compute effect-size estimates for all methods:
  estimates <- computeEstimates(persons)
  
  # Perform bootstrap and compute correlations:
  bootstrapEstimates <- lapply(1:settings$nBootstrap, performBootstrap, persons = persons)
  bootstrapEstimates <- bind_rows(bootstrapEstimates)
  correlations <- computeCorrelations(bootstrapEstimates)
  
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
    result <- tibble(logRr = logRr,
                     logCi95Lb = lb,
                     logCi95Ub = ub,
                     seLogRr = (ub - lb)/(2 * qnorm(0.975)))
    return(result)
  }
  averageModel <- computeAverage(estimates)
  averageModel$method <- "average"
  result <- bind_rows(estimates %>%
                        select("logRr", "logCi95Lb", "logCi95Ub", "seLogRr", "method") %>%
                        mutate(method = as.character(method)),
                      averageModel %>%
                        select("logRr", "logCi95Lb", "logCi95Ub", "seLogRr", "method"))
  result$iteration <- iteration
  return(result)
}

cluster <- makeCluster(10)
clusterRequire(cluster, "Cyclops")
clusterRequire(cluster, "dplyr")
clusterRequire(cluster, "tidyr")
results <- clusterApply(cluster, 1:100, runSimulation, settings = settings)
stopCluster(cluster)
results <- bind_rows(results)

# Compute coverage of the 95% CI per method (inc. model averaging):
for (method in unique(results$method)) {
  subset <- results[results$method == method, ]
  coverage <- mean(subset$logCi95Lb <= log(settings$rr) & subset$logCi95Ub >= log(settings$rr))
  writeLines(sprintf("Method: %s, coverage: %0.3f", 
                     method,
                     coverage))
}

