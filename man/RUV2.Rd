\name{RUV2}
\alias{RUV2}

\title{
Remove Unwanted Variation, 2-step
}

\description{
The RUV-2 algorithm.  Estimates and adjusts for unwanted variation using negative controls.
}

\usage{
RUV2(Y, X, ctl, k, Z=1, eta=NULL, include.intercept=TRUE,
fullW=NULL, svdyc=NULL, do_projectionplot=TRUE, inputcheck=TRUE)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{X}{
     The factor(s) of interest.  A m by p matrix, where m is the number of samples and p is the number of factors of interest.  Very often p = 1.  Factors and dataframes are also permissible, and converted to a matrix by \code{\link{design.matrix}}.
}
  \item{ctl}{
     An index vector to specify the negative controls.  Either a logical vector of length n or a vector of integers.  
}
  \item{k}{
     The number of unwanted factors to use.  Can be 0.
}
  \item{Z}{
     Any additional covariates to include in the model, typically a m by q matrix.  Factors and dataframes are also permissible, and converted to a matrix by \code{\link{design.matrix}}.  Alternatively, may simply be 1 (the default) for an intercept term.  May also be \code{NULL}.
}
  \item{eta}{
     Gene-wise (as opposed to sample-wise) covariates.  These covariates are adjusted for by RUV-1 before any further analysis proceeds.  Can be either (1) a matrix with n columns, (2) a matrix with n rows, (3) a dataframe with n rows, (4) a vector or factor of length n, or (5) simply 1, for an intercept term.  
}
  \item{include.intercept}{
     Applies to both \code{Z} and \code{eta}.  When \code{Z} or \code{eta} (or both) is specified (not \code{NULL}) but does not already include an intercept term, this will automatically include one.  If only one of \code{Z} or \code{eta} should include an intercept, this variable should be set to \code{FALSE}, and the intercept term should be included manually where desired.
}
  \item{fullW}{
     Can be included to speed up execution.  Is returned by previous calls of \code{RUV2} (see  below).
}
  \item{svdyc}{
     Can be included to speed up execution.  For internal use; please use \code{fullW} instead.
}
  \item{do_projectionplot}{
     Calculate the quantities necessary to generate a projection plot.
}
  \item{inputcheck}{
     Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
}
}

\details{
Implements the RUV-2 algorithm as described in Gagnon-Bartsch and Speed (2012), using the SVD as the factor analysis routine.  Unwanted factors W are estimated using control genes.  Y is then regressed on the variables X, Z, and W.   
}

\value{
A list containing
  \item{betahat}{The estimated coefficients of the factor(s) of interest.  A p by n matrix.}
  \item{sigma2}{Estimates of the features' variances.  A vector of length n.}
  \item{t}{t statistics for the factor(s) of interest.  A p by n matrix.}
  \item{p}{P-values for the factor(s) of interest.  A p by n matrix.}
  \item{Fstats}{F statistics for testing all of the factors in \code{X} simultaneously.}
  \item{Fpvals}{P-values for testing all of the factors in \code{X} simultaneously.}
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
  \item{eta}{\code{eta}.  Included for reference.}
  \item{fullW}{Can be used to speed up future calls of RUV2.}
  \item{projectionplotW}{A reparameterization of W useful for projection plots.}
  \item{projectionplotalpha}{A reparameterization of alpha useful for projection plots.}
  \item{include.intercept}{\code{include.intercept}.  Included for reference.}
  \item{method}{Character variable with value "RUV2".  Included for reference.}
}

\note{
Additional resources can be found at \href{http://www-personal.umich.edu/~johanngb/ruv/}{http://www-personal.umich.edu/~johanngb/ruv/}. 
}

\references{
Using control genes to correct for unwanted variation in microarray data.  Gagnon-Bartsch and Speed, 2012.  Available at: \href{http://biostatistics.oxfordjournals.org/content/13/3/539.full}{http://biostatistics.oxfordjournals.org/content/13/3/539.full}. 

Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
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

## Run RUV-2
fit = RUV2(Y, X, ctl, k)

## Get adjusted variances and p-values
fit = variance_adjust(fit)
}

\seealso{
\code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{variance_adjust}} 
}

\keyword{models}
\keyword{multivariate}



