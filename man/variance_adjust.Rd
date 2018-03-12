\name{variance_adjust}
\alias{variance_adjust}

\title{
Adjust Estimated Variances
}

\description{
Calculate rescaled variances, empirical variances, etc.  For use with RUV model fits.
}

\usage{
variance_adjust(fit, ctl.idx = NULL, ebayes = TRUE, pooled=TRUE, evar = TRUE,
                rsvar = TRUE, bin = 10, rescaleconst = NULL)
}

\arguments{
  \item{fit}{
  A RUV model fit (a list), as returned by RUV2 / RUV4 / RUVinv / RUVrinv
}
  \item{ctl.idx}{
  An index vector to specify the negative controls for use with the rescaled variances method.  If unspecified, by default \code{fit$ctl} is used.
}
  \item{ebayes}{
  A logical variable.  Should empirical bayes variance estimates be calculated?
}
  \item{pooled}{
  A logical variable.  Should pooled variance estimates be calculated?
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
    \item{sigma2.pooled}{Estimate of sigma^2 pooled (averaged) over all genes}
    \item{df.pooled}{Degrees of freedom for pooled estimate}
    \item{varbetahat}{"Standard" estimate of the variance of \code{betahat}}
    \item{varbetahat.rsvar}{"Rescaled Variances" estimate of the variance of \code{betahat}}
    \item{varbetahat.evar}{"Empirical Variances" estimate of the variance of \code{betahat}}
    \item{varbetahat.ebayes}{"Empirical Bayes" estimate of the variance of \code{betahat}}
    \item{varbetahat.rsvar.ebayes}{"Rescaled Empirical Bayes" estimate of the variance of \code{betahat}}
    \item{varbetahat.pooled}{"Pooled" estimate of the variance of \code{betahat}}
    \item{varbetahat.rsvar.pooled}{"Rescaled pooled" estimate of the variance of \code{betahat}}
    \item{varbetahat.evar.pooled}{Similar to the above, but all genes used to determine the rescaling, not just control genes}
    \item{p.rsvar}{P-values, after applying the method of rescaled variances}
    \item{p.evar}{P-values, after applying the method of empirical variances}
    \item{p.ebayes}{P-values, after applying the empirical bayes method of Smyth (2004)}
    \item{p.pooled}{P-values, after pooling variances}
    \item{p.rsvar.ebayes}{P-values, after applying the empirical bayes method of Smyth (2004) and the method of rescaled variances}
    \item{p.rsvar.pooled}{P-values, after pooling variances and the method of rescaled variances}
    \item{p.evar.pooled}{Similar to the above, but all genes used to determine the rescaling, not just control genes}
    \item{Fpvals.ebayes}{F test p-values, after applying the empirical bayes method of Smyth (2004)}
    \item{Fpvals.pooled}{F test p-values, after pooling variances}
    \item{p.BH}{FDR-adjusted p-values}
    \item{Fpvals.BH}{FDR-adjusted p-values (from F test)}
    \item{p.rsvar.BH}{FDR-adjusted p-values (from p.rsvar)}
    \item{p.evar.BH}{FDR-adjusted p-values (from p.evar)}
    \item{p.ebayes.BH}{FDR-adjusted p-values (from p.ebayes)}
    \item{p.rsvar.ebayes.BH}{FDR-adjusted p-values (from p.rsvar.ebayes)}
    \item{Fpvals.ebayes.BH}{FDR-adjusted F test p-values (from Fpvals.ebayes)}
    \item{p.pooled.BH}{FDR-adjusted p-values (from p.pooled)}
    \item{p.rsvar.pooled.BH}{FDR-adjusted p-values (from p.rsvar.pooled)}
    \item{p.evar.pooled.BH}{FDR-adjusted p-values (from p.evar.pooled)}
    \item{Fpvals.pooled.BH}{FDR-adjusted F test p-values (from Fpvals.pooled)}
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


