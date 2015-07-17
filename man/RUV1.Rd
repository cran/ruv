\name{RUV1}
\alias{RUV1}

\title{
RUV-1
}

\description{
The RUV-1 algorithm.  Generally used as a preprocessing step to RUV-2, RUV-4, RUV-inv, or RUV-rinv.
}

\usage{
RUV1(Y, eta, ctl)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{eta}{
     Gene-wise (as oposed to sample-wise) covariates.  A matrix with n columns.
}
  \item{ctl}{
     The negative controls.  A logical vector of length n.  
}
}

\details{
Implements the RUV-1 algorithm as described in Gagnon-Bartsch, Jacob, and Speed (2013).  Most often this algorithm is not used directly, but rather is called from RUV-2, RUV-4, RUV-inv, or RUV-rinv.
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
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}
}




