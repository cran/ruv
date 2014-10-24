\name{inputcheck1}
\alias{inputcheck1}

\title{
Input Check One
}

\description{
Performs a basic sanity check on the arguments passed to RUV-2, RUV-4, RUV-inv, and RUV-rinv.  
}

\usage{
inputcheck1(Y, X, Z, ctl)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{X}{
     The factor(s) of interest.  A m by p matrix, where m is the number of samples and p is the number of factors of interest.  Very often p = 1.
}
  \item{Z}{
     Any additional covariates to include in the model.  Either a m by q matrix of covariates, or simply 1 (the default) for an intercept term.
}
  \item{ctl}{
     The negative controls.  A logical vector of length n.  
}
}

\value{
Returns NULL.  The function is only called to check for problems in the arguments and to issue warnings if any problems are found.
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}
}



