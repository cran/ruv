\name{ruv_svdplot}
\alias{ruv_svdplot}

\title{
RUV SVD Plot
}

\description{
  A generalization of a PC (principal component) plot.
}

\usage{
ruv_svdplot(Y.data, Y.space = NULL, info = NULL, k = c(1, 2), Z = 1, left = TRUE)
}

\arguments{
  \item{Y.data}{
  The data matrix.  Rows are observations and columns are features (e.g. genes).
}
  \item{Y.space}{
  Either a data matrix of the same dimension as \code{Y.data}, or the SVD of such a matrix, as returned by the \code{svd} function.  The singular vectors of this matrix define the space in which \code{Y.data} will be plotted.  If NULL, \code{Y.data} itself is used.
}
  \item{info}{
  Additional data to be included in the ggplot, which can be used for setting aesthetics such as color.  Converted to a dataframe, which should have a number of rows equal to the number of rows of \code{Y.data} (if \code{left=TRUE}) or the number of columns of \code{Y.data} (if \code{left=FALSE}).
}
  \item{k}{
  A numeric vector of length 2.  The singular vectors to be plotted.  Typically integers, but fractional values can also be specified.  For example, a value of 2.5 corresponds to the linear combination (singular vector 2) + (singular vector 3), rescaled to have unit length.  Similarly, a value of 2.2 corresponds to the (rescaled) linear combination 8*(singular vector 2) + 2*(singular vector 3), and -2.2 corresponds to the (rescaled) linear combination 8*(singular vector 2) - 2*(singular vector 3).  Note that the vectors defined by 2.2 and -2.8 are orthogonal to each other, as are those defined by 2.3 and -2.7, etc.
}
  \item{Z}{
  Any variables to regress out of \code{Y.data} as a preprocessing step.  May simply be 1 (the default) for an intercept term, i.e. the columns of Y are mean centered.  May also be \code{NULL}.  Similarly for \code{Y.space}, unless \code{Y.space} is already an SVD.
}
  \item{left}{
  Plot the left singular vectors (if TRUE) or the right singular vectors (if FALSE).
}
}

\details{
  When \code{Y.space = NULL} and \code{Z = 1} and the values of \code{k} are integers, this is a standard PC plot.
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


