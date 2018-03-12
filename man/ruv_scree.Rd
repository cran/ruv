\name{ruv_scree}
\alias{ruv_scree}

\title{
RUV Scree Plot
}

\description{
A scree plot (on the log scale)
}

\usage{
ruv_scree(Y = NULL, Z = 1, Y.svd = NULL)
}

\arguments{
  \item{Y}{
  The data matrix.  Rows are observations and columns are features (e.g. genes).  If not specified, Y.svd must be specified instead (which is faster).
}
  \item{Z}{
  Any variables to regress out of Y as a preprocessing step.  May simply be 1 (the default) for an intercept term, i.e. the columns of Y are mean centered.  May also be \code{NULL}.
}
  \item{Y.svd}{
  The SVD of Y, as returned by the svd function.  
}
}

\details{
  Because 0 cannot be plotted on a log scale, if any singular values are equal to 0, they will be changed to the minimum non-zero singular value and plotted in red.  Exception: singular values that are 0 as a result of regressing out Z are simply not plotted.
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


