\name{invvar}
\alias{invvar}

\title{
Inverse Method Variances
}

\description{
  Estimate the features' variances using the inverse method.  This function is usually called from \code{\link{RUVinv}} and not normally intended for stand-alone use.
}

\usage{
invvar(Y, ctl, XZ = NULL, eta = NULL, lambda = NULL,
       invsvd = NULL)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{ctl}{
     The negative controls.  A logical vector of length n.  
}
  \item{XZ}{
     A m by (p + q) matrix containing both the factor(s) of interest (X) and known covariates (Z).
}
  \item{eta}{
     Gene-wise (as opposed to sample-wise) covariates.  These covariates are adjusted for by RUV-1 before any further analysis proceeds.  A matrix with n columns.
}
  \item{lambda}{
     Ridge parameter.  If specified, the ridged inverse method will be used.
}
  \item{invsvd}{
     Can be included to speed up execution.  Generally used when calling invvar many times with different values of lambda.
}
}

\value{
A list containing
  \item{sigma2}{Estimates of the features' variances.  A vector of length n.}
  \item{df}{The "effective degrees of freedom"}
  \item{invsvd}{Can be used to speed up future calls of invvar.}
}

\references{
Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}

\seealso{
\code{\link{RUVinv}}, \code{\link{RUVrinv}}
}


