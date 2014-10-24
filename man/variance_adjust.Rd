\name{variance_adjust}
\alias{variance_adjust}

\title{
Adjust Estimated Variances
}

\description{
Calculate rescaled variances, empirical variances, etc.  For use with RUV model fits.
}

\usage{
variance_adjust(fit, ebayes = TRUE, evar = TRUE, rsvar = TRUE, 
                bin = 10, rescaleconst = NULL)
}

\arguments{
  \item{fit}{
  A RUV model fit (a list), as returned by RUV2 / RUV4 / RUVinv / RUVrinv
}
  \item{ebayes}{
  A logical variable.  Should empirical bayes estimates be calculated?
}
  \item{evar}{
  A logical variable.  Should empirical variance estimates be calculated?
}
  \item{rsvar}{
  A logical variable.  Should rescaled variance estimates be calculated?
}
  \item{bin}{
  The bin size to use when calculating empirical variances.
}
  \item{rescaleconst}{
  Can be used to speed up execution.  See \code{\link{get_empirical_variances}}.
}
}

\value{
An RUV model fit (a list).  In addition to the elements of the list returned by RUV2 / RUV4 / RUVinv / RUVrinv, the list will now contain:
    \item{sigma2.ebayes}{Estimates of sigma^2 using the empirical bayes shrinkage method of Smyth (2004)}
    \item{df.ebayes}{Estimate of degrees of freedom using the empirical bayes shrinkage method of Smyth (2004)}
    \item{varbetahat}{"Standard" estimate of the variance of \code{betahat}}
    \item{varbetahat.rsvar}{"Rescaled Variances" estimate of the variance of \code{betahat}}
    \item{varbetahat.evar}{"Empirical Variances" estimate of the variance of \code{betahat}}
    \item{varbetahat.ebayes}{"Empirical Bayes" estimate of the variance of \code{betahat}}
    \item{varbetahat.rsvar.ebayes}{"Rescaled Empirical Bayes" estimate of the variance of \code{betahat}}
    \item{p.rsvar}{P-values, after applying the method of rescaled variances}
    \item{p.evar}{P-values, after applying the method of empirical variances}
    \item{p.ebayes}{P-values, after applying the empirical bayes method of Smyth (2004)}
    \item{p.rsvar.ebayes}{P-values, after applying the empirical bayes method of Smyth (2004) and the method of rescaled variances}
    \item{p.BH}{FDR-adjusted p-values}
    \item{p.rsvar.BH}{FDR-adjusted p-values, after applying the method of rescaled variances}
    \item{p.evar.BH}{FDR-adjusted p-values, after applying the method of empirical variances}
    \item{p.ebayes.BH}{FDR-adjusted p-values, after applying the empirical bayes method of Smyth (2004)}
    \item{p.rsvar.ebayes.BH}{FDR-adjusted p-values, after applying the empirical bayes method of Smyth (2004) and the method of rescaled variances}
}

\references{
Using control genes to correct for unwanted variation in microarray data.  Gagnon-Bartsch and Speed, 2012.  Available at: \href{http://biostatistics.oxfordjournals.org/content/13/3/539.full}{http://biostatistics.oxfordjournals.org/content/13/3/539.full}. 

Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.

Linear models and empirical bayes methods for assessing differential expression in microarray experiments.  Smyth, 2004.
}

\author{
Johann Gagnon-Bartsch
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

## Run RUV-inv
fit = RUVinv(Y, X, ctl)

## Get adjusted variances and p-values
fit = variance_adjust(fit)
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{get_empirical_variances}}, \code{\link{sigmashrink}}
}


