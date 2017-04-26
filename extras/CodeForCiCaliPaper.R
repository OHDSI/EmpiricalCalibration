# This is the code that is included as supplementary information for our CI calibration paper:


fitSystematicErrorModel <- function(logRr, seLogRr, trueLogRr) {
  if (any(is.infinite(seLogRr))) {
    warning("Estimate(s) with infinite standard error detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.infinite(seLogRr)]
    logRr <- logRr[!is.infinite(seLogRr)]
    seLogRr <- seLogRr[!is.infinite(seLogRr)]
  }
  if (any(is.infinite(logRr))) {
    warning("Estimate(s) with infinite logRr detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.infinite(logRr)]
    seLogRr <- seLogRr[!is.infinite(logRr)]
    logRr <- logRr[!is.infinite(logRr)]
  }
  if (any(is.na(seLogRr))) {
    warning("Estimate(s) with NA standard error detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.na(seLogRr)]
    logRr <- logRr[!is.na(seLogRr)]
    seLogRr <- seLogRr[!is.na(seLogRr)]
  }
  if (any(is.na(logRr))) {
    warning("Estimate(s) with NA logRr detected. Removing before fitting error model")
    trueLogRr <- trueLogRr[!is.na(logRr)]
    seLogRr <- seLogRr[!is.na(logRr)]
    logRr <- logRr[!is.na(logRr)]
  }
  
  gaussianProduct <- function(mu1, mu2, sd1, sd2) {
    (2 * pi)^(-1/2) * (sd1^2 + sd2^2)^(-1/2) * exp(-(mu1 - mu2)^2/(2 * (sd1^2 + sd2^2)))
  }
  
  LL <- function(theta, logRr, seLogRr, trueLogRr) {
    result <- 0
    for (i in 1:length(logRr)) {
      mean <- theta[1] + theta[2] * trueLogRr[i]
      sd <- exp(theta[3] + theta[4] * trueLogRr[i])
      result <- result - log(gaussianProduct(logRr[i], mean, seLogRr[i], sd))
    }
    if (is.infinite(result))
      result <- 99999
    result
  }
  theta <- c(0, 1, -2, 0)
  fit <- optim(theta,
               LL,
               logRr = logRr,
               seLogRr = seLogRr,
               trueLogRr = trueLogRr,
               method = "BFGS",
               hessian = TRUE,
               control = list(parscale = c(1, 1, 10, 10)))
  fisher_info <- solve(fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  model <- fit$par
  names(model) <- c("meanIntercept", "meanSlope", "logSdIntercept", "logSdSlope")
  fisher_info <- solve(fit$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  attr(model, "CovarianceMatrix") <- fisher_info
  attr(model, "LB95CI") <- fit$par + qnorm(0.025) * prop_sigma
  attr(model, "UB95CI") <- fit$par + qnorm(0.975) * prop_sigma
  class(model) <- "systematicErrorModel"
  model
}

calibrateConfidenceInterval <- function(logRr, seLogRr, model, ciWidth = 0.95) {
  
  opt <- function(x,
                  z,
                  logRr,
                  se,
                  interceptMean,
                  slopeMean,
                  interceptLogSd,
                  slopeLogSd) {
    mean <- interceptMean + slopeMean * x
    sd <- exp(interceptLogSd + slopeLogSd * x)
    return(z + (mean - logRr)/sqrt((sd)^2 + (se)^2))
  }
  
  logBound <- function(ciWidth,
                       lb = TRUE,
                       logRr,
                       se,
                       interceptMean,
                       slopeMean,
                       interceptLogSd,
                       slopeLogSd) {
    z <- qnorm((1 - ciWidth)/2)
    if (lb) {
      z <- -z
    }
    # Simple grid search for upper bound where opt is still positive:
    upper <- -9
    while(opt(x = upper,
              z = z,
              logRr = logRr,
              se = se,
              interceptMean = interceptMean,
              slopeMean = slopeMean,
              interceptLogSd = interceptLogSd,
              slopeLogSd = slopeLogSd) < 0) {
      upper <- upper + 1
    }
    
    uniroot(f = opt,
            interval = c(-10, upper),
            z = z,
            logRr = logRr,
            se = se,
            interceptMean = interceptMean,
            slopeMean = slopeMean,
            interceptLogSd = interceptLogSd,
            slopeLogSd = slopeLogSd)$root
  }
  
  result <- data.frame(logRr = rep(0, length(logRr)), logLb95Rr = 0, logUb95Rr = 0)
  for (i in 1:nrow(result)) {
    if (is.infinite(logRr[i]) || is.na(logRr[i]) || is.infinite(seLogRr[i]) || is.na(seLogRr[i])) {
      result$logRr[i] <- NA
      result$logLb95Rr[i] <- NA
      result$logUb95Rr[i] <- NA
    } else {
      result$logRr[i] <- logBound(0,
                                  TRUE,
                                  logRr[i],
                                  seLogRr[i],
                                  model[1],
                                  model[2],
                                  model[3],
                                  model[4])
      result$logLb95Rr[i] <- logBound(ciWidth,
                                      TRUE,
                                      logRr[i],
                                      seLogRr[i],
                                      model[1],
                                      model[2],
                                      model[3],
                                      model[4])
      result$logUb95Rr[i] <- logBound(ciWidth,
                                      FALSE,
                                      logRr[i],
                                      seLogRr[i],
                                      model[1],
                                      model[2],
                                      model[3],
                                      model[4])
    }
  }
  result$seLogRr <- (result$logLb95Rr - result$logUb95Rr)/(2 * qnorm((1 - ciWidth)/2))
  return(result)
}