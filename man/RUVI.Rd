\name{RUVI}
\alias{RUVI}
\alias{RUV1}

\title{
RUV-I
}

\description{
The RUV-I algorithm.  Generally used as a preprocessing step to RUV-2, RUV-4, RUV-inv, RUV-rinv, or RUVIII.  RUV1 is an alias of (identical to) RUVI.
}

\usage{
RUVI(Y, eta, ctl, include.intercept = TRUE)

RUV1(Y, eta, ctl, include.intercept = TRUE)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{eta}{
     Gene-wise (as opposed to sample-wise) covariates.  A matrix with n columns.
}
  \item{ctl}{
     The negative controls.  A logical vector of length n.  
}
  \item{include.intercept}{
     Add an intercept term to eta if it does not include one already.
}
}

\details{
Implements the RUV-I algorithm as described in Gagnon-Bartsch, Jacob, and Speed (2013).  Most often this algorithm is not used directly, but rather is called from RUV-2, RUV-4, RUV-inv, or RUV-rinv.  Note that RUV1 and RUVI are two different names for the same (identical) function.
}

\value{
An adjusted data matrix (i.e., an adjusted Y)
}

\references{
Using control genes to correct for unwanted variation in microarray data.  Gagnon-Bartsch and Speed, 2012.  Available at: \href{http://biostatistics.oxfordjournals.org/content/13/3/539.full}{http://biostatistics.oxfordjournals.org/content/13/3/539.full}. 

Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{RUVIII}}
}




