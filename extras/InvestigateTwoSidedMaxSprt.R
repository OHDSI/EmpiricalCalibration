# Scenario: we want to monitor an exposure for an outcome over time, as data accrues.
# We want to keep on testing for an increased risk, until we either conclude there's
# an effect, or we concluded it is safe.

# Simple Poisson process, with known expected rate ----------------------------------------------
library(EmpiricalCalibration)

alpha <- 0.05
beta <- 0.20
theta1 <- log(1.25) # Minimally clinically relevant increased risk


# How much extra sample (in expected count) at each look?
getIncrementalExpectedCount <- function(t) {
  return(10)
}

computeLlr <- function(observed, expected, theta1) {
  rr <- exp(theta1)
  if (observed < expected) {
    ll0 <-  dpois(observed, observed, log = TRUE)
  } else {
    ll0 <-  dpois(observed, expected, log = TRUE)
  }
  if (observed > expected * rr) {
    ll1 <-  dpois(observed, observed, log = TRUE)
  } else {
    ll1 <-  dpois(observed, expected * rr, log = TRUE)
  }
  llr <- ll1 - ll0
  return(llr)
}

# Using crude grid + Monte-Carlo for now:
computeCvs <- function(alpha, beta, theta1, getIncrementalExpectedCount, nSimulations = 1000, maxT = 1000) {

  evaluateCvs <- function(cvLower, cvUpper, irr) {
    signals <- 0
    safe <- 0
    for (i in 1:nSimulations) {
      observed <- 0
      expected <- 0
      for (t in 1:maxT) {
        groupSize <- getIncrementalExpectedCount(t)
        expected <- expected + groupSize
        observed <- observed + rpois(1, groupSize * irr)
        llr <- computeLlr(observed, expected, theta1)
        if (llr > cvUpper) {
          signals <- signals + 1
          break
        }
        if (llr < cvLower) {
          safe <- safe + 1
          break
        }
      }
    }
    return(c(positiveRate = signals / nSimulations,
             negativeRate = safe / nSimulations))
  }
  # From Wald, 1945:
  cvUpperCenter <- log((1 - beta) / alpha)
  cvLowerCenter <- log(beta / (1 - alpha))
  
  cvLower <- cvLowerCenter + seq(from = -0.5, to = 1.5, by = 0.1)
  cvUpper <- cvUpperCenter + seq(from = -1.5, to = 0.5, by = 0.1)
  grid <- expand.grid(cvLower = cvLower, cvUpper = cvUpper)
  for (i in 1:nrow(grid)) {
    print(sprintf("Evaluating %d of %d", i, nrow(grid)))
    grid$type1[i] <- evaluateCvs(cvLower = grid$cvLower[i], cvUpper = grid$cvUpper[i], irr = 1)[1]
    grid$type2[i] <- evaluateCvs(cvLower = grid$cvLower[i], cvUpper = grid$cvUpper[i], irr = exp(theta1))[2]
  }
  grid$distance <- sqrt((grid$type1 - alpha) ^ 2 + (grid$type2 - beta) ^ 2)
  
  # library(ggplot2)
  # ggplot(grid, aes(x = cvLower, y = cvUpper, fill = distance)) +
  #   geom_tile()
  
  winner <- which(grid$distance == min(grid$distance))[1]
  return(grid[winner, ])
}
# Compute critical values:
cvs <- computeCvs(alpha, beta, theta1, getIncrementalExpectedCount)

# Try out critical values:
cvUpper <- cvs$cvUpper
cvLower <- cvs$cvLower
trueRr <- 1
nSimulations <- 1000
maxT <- 1000

signals <- 0
safe <- 0
timeToDone <- rep(0, nSimulations)
for (i in 1:nSimulations) {
  observed <- 0
  expected <- 0
  for (t in 1:maxT) {
    groupSize <- getIncrementalExpectedCount(t)
    expected <- expected + groupSize
    observed <- observed + rpois(1, groupSize * trueRr)
    llr <- computeLlr(observed, expected, theta1)
    if (llr > cvUpper) {
      signals <- signals + 1
      break
    }
    if (llr < cvLower) {
      safe <- safe + 1
      break
    }
  }
  timeToDone[i] <- t
}
writeLines(sprintf("True RR: %0.2f, Positives: %0.1f%%, Negatives: %0.1f%%, Min time: %d, Median time: %0.1f, Max time: %d", 
                   trueRr,
                   100 * signals / nSimulations,
                   100 * safe / nSimulations,
                   min(timeToDone),
                   median(timeToDone),
                   max(timeToDone)))

