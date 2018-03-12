\name{projectionplotvariables}
\alias{projectionplotvariables}

\title{
Projection Plot Variables
}

\description{
Calculates the variables necessary to produce a projection plot.
}

\usage{
projectionplotvariables(Y, X, W)
}

\arguments{
  \item{Y}{
    The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{X}{
     A m by p matrix containing the factor(s) of interest.
}
  \item{W}{
     A m by k matrix containing the estimated unwanted factors.
}
}

\details{
Typically intended for internal use, and called after adjustment for known covariates (Z).
}

\value{
A list containing
  \item{byx}{Regression coefficients from regressing Y on X.}
  \item{bwx}{Regression coefficients from regressing W on X.}
  \item{projectionplotalpha}{A reparameterization of alpha.}
  \item{projectionplotW}{A reparameterization of W.}
}


