\name{RUV4}
\alias{RUV4}

\title{
Remove Unwanted Variation, 4-step
}

\description{
The RUV-4 algorithm.  Estimates and adjusts for unwanted variation using negative controls.
}

\usage{
RUV4(Y, X, ctl, k, Z = 1, eta = NULL, fullW0 = NULL,
     inputcheck = TRUE)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{X}{
     The factor(s) of interest.  A m by p matrix, where m is the number of samples and p is the number of factors of interest.  Very often p = 1.
}
  \item{ctl}{
     The negative controls.  A logical vector of length n.  
}
  \item{k}{
     The number of unwanted factors to use.  Can be 0.
}
  \item{Z}{
     Any additional covariates to include in the model.  Either a m by q matrix of covariates, or simply 1 (the default) for an intercept term.
}
  \item{eta}{
     Gene-wise (as oposed to sample-wise) covariates.  These covariates are adjusted for by RUV-1 before any further analysis proceeds.  A matrix with n columns.
}
  \item{fullW0}{
     Can be included to speed up execution.
}
  \item{inputcheck}{
     Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
}
}

\details{
Implements the RUV-4 algorithm as described in Gagnon-Bartsch, Jacob, and Speed (2013), using the SVD as the factor analysis routine.  Unwanted factors W are estimated using control genes.  Y is then regressed on the variables X, Z, and W.   
}

\value{
A list containing
  \item{betahat}{The estimated coefficients of the factor(s) of interest.  A p by n matrix.}
  \item{sigma2}{Estimates of the features' variances.  A vector of length n.}
  \item{t}{t statistics for the factor(s) of interest.  A p by n matrix.}
  \item{p}{P-values for the factor(s) of interest.  A p by n matrix.}
  \item{multiplier}{The constant by which \code{sigma2} must be multiplied in order get an estimate of the variance of \code{betahat}}
  \item{df}{The number of residual degrees of freedom.}
  \item{W}{The estimated unwanted factors.}
  \item{alpha}{The estimated coefficients of W.}
  \item{byx}{The coefficients in a regression of Y on X (after both Y and X have been "adjusted" for Z).  Useful for projection plots.}
  \item{bwx}{The coefficients in a regression of W on X (after X has been "adjusted" for Z).  Useful for projection plots.}
  \item{X}{\code{X}.  Included for reference.}
  \item{k}{\code{k}.  Included for reference.}
  \item{ctl}{\code{ctl}.  Included for reference.}
  \item{Z}{\code{Z}.  Included for reference.}
  \item{fullW0}{Can be used to speed up future calls of RUV4.}
}

\note{
Additional resources can be found at \href{http://www.stat.berkeley.edu/~johann/ruv}{http://www.stat.berkeley.edu/~johann/ruv}. 
}

\references{
Using control genes to correct for unwanted variation in microarray data.  Gagnon-Bartsch and Speed, 2012.  Available at: \href{http://biostatistics.oxfordjournals.org/content/13/3/539.full}{http://biostatistics.oxfordjournals.org/content/13/3/539.full}. 

Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.
}

\author{
Johann Gagnon-Bartsch \email{johann@stat.berkeley.edu}
}

\examples{
## Create some simulated data
m = 50
n = 10000
nc = 1000
p = 1
k = 20
ctl = rep(FALSE, n)
ctl[1:nc] = TRUE
X = matrix(c(rep(0,floor(m/2)), rep(1,ceiling(m/2))), m, p)
beta = matrix(rnorm(p*n), p, n)
beta[,ctl] = 0
W = matrix(rnorm(m*k),m,k)
alpha = matrix(rnorm(k*n),k,n)
epsilon = matrix(rnorm(m*n),m,n)
Y = X\%*\%beta + W\%*\%alpha + epsilon

## Run RUV-4
fit = RUV4(Y, X, ctl, k)

## Get adjusted variances and p-values
fit = variance_adjust(fit)
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{variance_adjust}} 
}

\keyword{models}
\keyword{multivariate}






