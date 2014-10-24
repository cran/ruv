\name{RUVrinv}
\alias{RUVrinv}

\title{
Remove Unwanted Variation, ridged inverse method
}

\description{
The RUV-rinv algorithm.  Estimates and adjusts for unwanted variation using negative controls.
}

\usage{
RUVrinv(Y, X, ctl, Z = 1, eta = NULL, fullW0 = NULL,
        invsvd = NULL, lambda = NULL, k = NULL, 
        l = NULL, randomization = FALSE, iterN = 1e+05,
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
  \item{Z}{
     Any additional covariates to include in the model.  Either a m by q matrix of covariates, or simply 1 (the default) for an intercept term.
}
  \item{eta}{
     Gene-wise (as oposed to sample-wise) covariates.  These covariates are adjusted for by RUV-1 before any further analysis proceeds.  A matrix with n columns.
}
  \item{fullW0}{
     Can be included to speed up execution.
}
  \item{invsvd}{
     Can be included to speed up execution.  Generally used when calling RUV(r)inv many times with different values of lambda.
}
  \item{lambda}{
     Ridge parameter.  If unspecified, an appropriate default will be used.
}
  \item{k}{
     When calculating the default value of lambda, a call to RUV4 is made.  This parameter specifies the value of k to use.  Otherwise, an appropriate default k will be used.
}
  \item{l}{
     If lambda and k are both NULL, then k must be estimated using the getK routine.  The getK routine only accepts a single-column X.  If p > 1, l specifies which column of X should be used in the getK routine. 
}
  \item{randomization}{
     Whether the inverse-method variances should be computed using randomly generated factors of interest (as opposed to a numerical integral).
}
  \item{iterN}{
     The number of random "factors of interest" to generate (used only when randomization=TRUE).
}
  \item{inputcheck}{
     Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
}
}

\details{
Implements the RUV-rinv algorithm as described in Gagnon-Bartsch, Jacob, and Speed (2013).  This function is essentially just a wrapper to RUVinv, but with a little extra code to calculate the default value of \code{lambda}.  
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
  \item{lambda}{\code{lambda}.  Included for reference.}
  \item{invsvd}{Can be used to speed up future calls of RUV(r)inv.}
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

## Run RUV-rinv
fit = RUVrinv(Y, X, ctl)

## Get adjusted variances and p-values
fit = variance_adjust(fit)
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{variance_adjust}}, \code{\link{invvar}}, \code{\link{getK}}  
}

\keyword{models}
\keyword{multivariate}


