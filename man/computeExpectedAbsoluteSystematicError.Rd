% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExpectedSystematicError.R
\name{computeExpectedAbsoluteSystematicError}
\alias{computeExpectedAbsoluteSystematicError}
\title{Compute the expected absolute systematic error}
\usage{
computeExpectedAbsoluteSystematicError(null, alpha = 0.05)
}
\arguments{
\item{null}{An object of class \code{null} created using the \code{fitNull} function or an
object of class \code{mcmcNull} created using the \code{fitMcmcNull} function.}

\item{alpha}{The expected type I error for computing the credible interval.}
}
\value{
The expected absolute systematic error. If the provided \code{null} argument is of type \code{mcmcNull},
the credible interval (defined by \code{alpha}) is also returned.
}
\description{
For a random study estimate, what is the expected value of the absolute systematic error?
Provides a single summary value for a null distribution. The expected systematic error of a null
distribution is equal to its mean (mu), and is insensitive to the spread of the null distribution (sigma).

Taking the absolute value of the expected systematic error we can express both mean and spread of the
estimated null distribution.
}
\examples{
data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]
null <- fitNull(negatives$logRr, negatives$seLogRr)
computeExpectedAbsoluteSystematicError(null)
}
\seealso{
\code{\link{compareEase}} for comparing the expected absolute systematic error of two sets of estimates for the same negative controls.
}
