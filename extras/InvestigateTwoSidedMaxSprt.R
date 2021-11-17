# Scenario: we want to monitor an exposure for an outcome over time, as data accrues.
# We want to keep on testing for an increased risk, until we either conclude there's
# an effect, or we concluded it is safe.

# Simple Poisson process, with known expected rate ----------------------------------------------
library(EmpiricalCalibration)

alpha <- 0.05
beta <- 0.20
rr1 <- 1.25 # Minimally clinically relevant increased risk


# How much extra sample (in expected count) at each look?
getIncrementalExpectedCount <- function(t) {
  return(10)
}

computeLlr <- function(observed, expected, rr1) {
  if (observed < expected) {
    ll0 <-  dpois(observed, observed, log = TRUE)
  } else {
    ll0 <-  dpois(observed, expected, log = TRUE)
  }
  if (observed > expected * rr1) {
    ll1 <-  dpois(observed, observed, log = TRUE)
  } else {
    ll1 <-  dpois(observed, expected * rr, log = TRUE)
  }
  llr <- ll1 - ll0
  return(llr)
}

# Using crude grid + Monte-Carlo for now:
computeCvs <- function(alpha, beta, rr1, getIncrementalExpectedCount, nSimulations = 1000, maxT = 1000) {

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
        llr <- computeLlr(observed, expected, rr1)
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
    message(sprintf("Evaluating %d of %d", i, nrow(grid)))
    grid$type1[i] <- evaluateCvs(cvLower = grid$cvLower[i], cvUpper = grid$cvUpper[i], irr = 1)[1]
    grid$type2[i] <- evaluateCvs(cvLower = grid$cvLower[i], cvUpper = grid$cvUpper[i], irr = rr1)[2]
  }
  grid$distance <- sqrt((grid$type1 - alpha) ^ 2 + (grid$type2 - beta) ^ 2)
  
  # library(ggplot2)
  # ggplot(grid, aes(x = cvLower, y = cvUpper, fill = distance)) +
  #   geom_tile()
  
  winner <- which(grid$distance == min(grid$distance))[1]
  return(grid[winner, ])
}
# Compute critical values:
cvs <- computeCvs(alpha, beta, rr1, getIncrementalExpectedCount)

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
    llr <- computeLlr(observed, expected, rr1)
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



# Binomial process ----------------------------------------------------------------
library(EmpiricalCalibration)

alpha <- 0.05
beta <- 0.20
rr1 <- 1.5 # Minimally clinically relevant increased risk
z <- 1 # Exposed / unexposed

# How much extra sample (in expected count) at each look?
getIncrementalExpectedCount <- function(t) {
  return(10)
}

computeLlr <- function(observed, count, rr1, z) {
  pNull <- 1 / (1 + z)
  pAlternative <- rr1 / (rr1 + z)
  observedP <- observed / count
  if (observedP < pNull) {
    ll0 <- dbinom(observed, count, observedP, log = TRUE)
  } else {
    ll0 <- dbinom(observed, count, pNull, log = TRUE)
  }
  if (observedP > pAlternative) {
    ll1 <- dbinom(observed, count, observedP, log = TRUE)
  } else {
    ll1 <- dbinom(observed, count, pAlternative, log = TRUE)
  }
  llr <- ll1 - ll0
  return(llr)
}

# Using crude grid + Monte-Carlo for now:
computeCvs <- function(alpha, beta, rr1, getIncrementalExpectedCount, nSimulations = 1000, maxT = 1000) {
  
  evaluateCvs <- function(cvLower, cvUpper, irr) {
    signals <- 0
    safe <- 0
    for (i in 1:nSimulations) {
      count <- 0
      observed <- 0
      for (t in 1:maxT) {
        groupSize <- getIncrementalExpectedCount(t)
        count <- count + groupSize
        observed <- observed + rbinom(1, groupSize, irr / (irr + z))
        llr <- computeLlr(observed, count, rr1, z)
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
    message(sprintf("Evaluating %d of %d", i, nrow(grid)))
    grid$type1[i] <- evaluateCvs(cvLower = grid$cvLower[i], cvUpper = grid$cvUpper[i], irr = 1)[1]
    grid$type2[i] <- evaluateCvs(cvLower = grid$cvLower[i], cvUpper = grid$cvUpper[i], irr = rr1)[2]
  }
  grid$distance <- sqrt((grid$type1 - alpha) ^ 2 + (grid$type2 - beta) ^ 2)
  
  # library(ggplot2)
  # ggplot(grid, aes(x = cvLower, y = cvUpper, fill = distance)) +
  #   geom_tile()
  
  winner <- which(grid$distance == min(grid$distance))[1]
  return(grid[winner, ])
}
# Compute critical values:
cvs <- computeCvs(alpha, beta, rr1, getIncrementalExpectedCount)

# Try out critical values:
cvLower = -1.658145; cvUpper = 3.272589

cvUpper <- cvs$cvUpper
cvLower <- cvs$cvLower
trueRr <- 4
nSimulations <- 1000
maxT <- 1000

signals <- 0
safe <- 0
timeToDone <- rep(0, nSimulations)
for (i in 1:nSimulations) {
  observed <- 0
  count <- 0
  for (t in 1:maxT) {
    groupSize <- getIncrementalExpectedCount(t)
    count <- count + groupSize
    observed <- observed + rbinom(1, groupSize, trueRr / (trueRr + z))
    llr <- computeLlr(observed, count, rr1, z)
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
# rr1 = 1.25
# True RR: 1.00, Positives: 5.5%, Negatives: 94.5%, Min time: 1, Median time: 13.0, Max time: 205
# True RR: 1.25, Positives: 77.3%, Negatives: 22.7%, Min time: 1, Median time: 23.5, Max time: 240
# True RR: 1.50, Positives: 96.8%, Negatives: 3.2%, Min time: 1, Median time: 12.0, Max time: 65
# True RR: 2.00, Positives: 97.6%, Negatives: 2.4%, Min time: 1, Median time: 5.0, Max time: 25
# True RR: 4.00, Positives: 99.9%, Negatives: 0.1%, Min time: 1, Median time: 2.0, Max time: 7
# 
# rr1 = 1.5
# True RR: 1.00, Positives: 4.8%, Negatives: 95.2%, Min time: 1, Median time: 5.0, Max time: 61
# True RR: 1.25, Positives: 43.4%, Negatives: 56.6%, Min time: 1, Median time: 10.0, Max time: 73
# True RR: 1.50, Positives: 83.5%, Negatives: 16.5%, Min time: 1, Median time: 7.0, Max time: 83
# True RR: 2.00, Positives: 96.2%, Negatives: 3.8%, Min time: 1, Median time: 5.0, Max time: 21
# True RR: 4.00, Positives: 100.0%, Negatives: 0.0%, Min time: 1, Median time: 2.0, Max time: 8
