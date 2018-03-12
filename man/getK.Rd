\name{getK}
\alias{getK}

\title{
  Get K
}

\description{
  Finds an often-suitable value of K for use in RUV-4.
}

\usage{
getK(Y, X, ctl, Z = 1, eta = NULL, include.intercept = TRUE,
     fullW0 = NULL, cutoff = NULL, method="select", l=1, inputcheck = TRUE)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of samples and n is the number of features.
}
  \item{X}{
     The factor(s) of interest.  A m by p matrix, where m is the number of samples and p is the number of factors of interest.  Note that X should be only a single column, i.e. p = 1; if X has more than one column, only column \code{l} will be used (see below).  
}
  \item{ctl}{
     An index vector to specify the negative controls.  Either a logical vector of length n or a vector of integers.  
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
  \item{fullW0}{
     Can be included to speed up execution.  Is returned by previous calls of \code{getK}, \code{RUV4}, \code{RUVinv}, or \code{RUVrinv} (see  below).
}
  \item{cutoff}{
     Specify an alternative cut-off.  Default is the (approximate) 90th percentile of the distribution of the first singular value of an m by n gaussian matrix.
}
  \item{method}{
     Can be set to either \code{leave1out}, \code{fast}, or \code{select}.  \code{leave1out} is the preferred method but may be slow, \code{fast} is an approximate method that is faster but may provide poor results if n_c is not much larger than m, and \code{select} (the default) tries to choose for you.
}
  \item{l}{
     Which column of X to use in the getK algorithm. 
}
  \item{inputcheck}{
     Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
}
}

\value{
A list containing
  \item{k}{the estimated value of k}
  \item{cutoff}{The cutoff value used}
  \item{sizeratios}{A measure of the relative sizes of the rows of alpha.}
  \item{fullW0}{Can be used to speed up future calls of RUV4.}
}

\section{Warning }{This value of K will not be the best choice in all cases.  Moreover, it will often be a poor choice of K for use with RUV2.  See Gagnon-Bartsch and Speed (2012) for commentary on estimating k.}

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

## Run getK
temp = getK(Y, X, ctl)
K = temp$k
}


\seealso{
\code{\link{RUV4}}
}



