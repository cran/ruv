\name{RUVIII}
\alias{RUVIII}

\title{
RUV-III
}

\description{
Globally adjust data matrix using both negative controls and replicates.
}

\usage{
RUVIII(Y, M, ctl, k = NULL, eta = NULL, include.intercept = TRUE,
       average = FALSE, fullalpha = NULL, return.info = FALSE, inputcheck = TRUE)
}

\arguments{
  \item{Y}{
     The data.  A m by n matrix, where m is the number of observations and n is the number of features.
}
  \item{M}{
     The replicate structure.  Represented internally as a mapping matrix.  The mapping matrix has m rows (one for each observation), and each column represents a set of replicates.  The (i, j)-th entry of the mapping matrix is 1 if the i-th observation is in replicate set j, and 0 otherwise.  Each observation must be in exactly one set of replicates (some replicate sets may contain only one observation), and thus each row of M must sum to 1.  \code{M} may be the mapping matrix itself.  Alternatively, \code{M} may be a vector, factor, or dataframe, in which case it is converted to the mapping matrix by the \code{\link{replicate.matrix}} function.
}
  \item{ctl}{
     An index vector to specify the negative controls.  Either a logical vector of length n or a vector of integers.
}
  \item{k}{
     The number of unwanted factors to use.  Can be 0, in which case no adjustment is made.  Can also be NULL (the default value), in which case the maximum possible value of k is used; note that in this case no singular value decomposition is necessary and execution is faster.
}
  \item{eta}{
     Gene-wise (as opposed to sample-wise) covariates.  These covariates are adjusted for by RUV-1 before any further analysis proceeds.  Can be either (1) a matrix with n columns, (2) a matrix with n rows, (3) a dataframe with n rows, (4) a vector or factor of length n, or (5) simply 1, for an intercept term.  
}
  \item{include.intercept}{
     When \code{eta} is specified (not \code{NULL}) but does not already include an intercept term, this will automatically include one.  
}
  \item{average}{
     Average replicates after adjustment. 
}
  \item{fullalpha}{
     Can be included to speed up execution.  Is returned by previous calls of \code{RUVIII} (see  below).
}
  \item{return.info}{
     If \code{FALSE}, only the adjusted data matrix is returned.  If \code{TRUE}, additional information is returned (see below).
}
  \item{inputcheck}{
     Perform a basic sanity check on the inputs, and issue a warning if there is a problem.
}
}

\value{
If code{return.info} is \code{TRUE}, a list is returned that contains:
  \item{newY}{The adjusted data matrix.}
  \item{M}{The replicate mapping matrix.  Included for reference.}
  \item{fullalpha}{Can be used to speed up future calls to \code{RUVIII}}
Otherwise, if \code{return.info} is \code{FALSE}, only the adjusted data matrix is returned.
}

\note{
Additional resources can be found at \href{http://www-personal.umich.edu/~johanngb/ruv/}{http://www-personal.umich.edu/~johanngb/ruv/}. 
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}

\keyword{models}
\keyword{multivariate}
