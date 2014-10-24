\name{get_empirical_variances}
\alias{get_empirical_variances}

\title{
Get empirical variances
}

\description{
This method implements the method of empirical variances as described in Gagnon-Bartsch, Jacob, and Speed (2013).  This function is normally called from the function \code{\link{variance_adjust}}, and is not normally intended for stand-alone use.
}

\usage{
get_empirical_variances(sigma2, betahat, bin = 10,
                        rescaleconst = NULL)
}

\arguments{
  \item{sigma2}{
  Estimates of sigma^2
}
  \item{betahat}{
  Estimates of beta
}
  \item{bin}{
  The bin size
}
  \item{rescaleconst}{
  The expected value of the average of the smallest \code{bin - 1} of \code{bin} independent chi-square random variables.  This can be specified to save computational time (otherwise, it is calculated by simulation).
}
}


\value{
A vector of the empirical variances.
}

\references{
Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.
}

\author{
Johann Gagnon-Bartsch
}

\seealso{
\code{\link{variance_adjust}}
}



