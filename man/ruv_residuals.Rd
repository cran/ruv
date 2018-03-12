\name{ruv_residuals}
\alias{ruv_residuals}

\title{
RUV Residuals
}

\description{
Calculate the residuals or adjusted data matrix of an RUV2 or RUV4 fit.
}

\usage{
ruv_residuals(fit, type=c("residuals", "adjusted.Y"), subset_and_sort=TRUE)
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{type}{
  Whether to compute residuals or an adjusted data matrix.  Caution; see details below.
}
  \item{subset_and_sort}{
  Whether to subset and sort the features, as in \code{\link{ruv_summary}}.
}
}

\details{
  This function will return either the residuals or an adjusted data matrix.  The residuals are the result of removing all factors (wanted and unwanted), whereas the adjusted data matrix is the result of removing only the unwanted factors.  
  
  The residuals can be useful for diagnostics, e.g. in producing a residual RLE plot.  The adjusted data matrix may also be useful for diagnostics, but typically should *not* be used for any additional downstream analyses.  The adjusted data matrix can suffer from overfitting, which can be severe, especially when k is large, and this can produce artificially "good" results in downstream analyses.  
  
  If an adjusted data matrix for use in downstream analyses is desired, see \code{\link{RUVIII}}.
}

\value{
Either a matrix of residuals, or an adjusted data matrix.
}

\author{
Johann Gagnon-Bartsch
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{ruv_summary}}, \code{\link{RUVIII}}
}

