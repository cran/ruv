\name{ruv_svdgridplot}
\alias{ruv_svdgridplot}

\title{
RUV SVD Grid Plot
}

\description{
  A plot composed of a grid of several subplots created by \code{\link{ruv_svdplot}}
}

\usage{
ruv_svdgridplot(Y.data, Y.space = NULL, rowinfo = NULL, colinfo = NULL, k = 1:3, Z = 1,
                left.additions = NULL, right.additions = NULL, 
                factor.labels = paste("S.V.", k))
}

\arguments{
  \item{Y.data}{
  The data matrix.  Rows are observations and columns are features (e.g. genes).
}
  \item{Y.space}{
  Either a data matrix of the same dimension as \code{Y.data}, or the SVD of such a matrix, as returned by the \code{svd} function.  The singular vectors of this matrix define the space in which \code{Y.data} will be plotted.  If NULL, \code{Y.data} itself is used.
}
  \item{rowinfo}{
  A dataframe of information about the observations.  Should have the same number of rows as Y.  This information will be included in the ggplots, and can be used for setting aesthetics such as color.
}
  \item{colinfo}{
  A dataframe of information about the observations.  Should have a number of rows equal to the number of columns of Y.  This information will be included in the ggplots, and can be used for setting aesthetics such as color.
}
  \item{k}{
  A numeric vector of the singular vectors to be plotted.  Typically integers, but fractional values can also be specified.  For example, a value of 2.5 corresponds to the linear combination (singular vector 2) + (singular vector 3), rescaled to have unit length.  Similarly, a value of 2.2 corresponds to the (rescaled) linear combination 8*(singular vector 2) + 2*(singular vector 3), and -2.2 corresponds to the (rescaled) linear combination 8*(singular vector 2) - 2*(singular vector 3).  Note that the vectors defined by 2.2 and -2.8 are orthogonal to each other, as are those defined by 2.3 and -2.7, etc.
}
  \item{Z}{
  Any variables to regress out of \code{Y.data} as a preprocessing step.  May simply be 1 (the default) for an intercept term, i.e. the columns of Y are mean centered.  May also be \code{NULL}.  Similarly for \code{Y.space}, unless \code{Y.space} is already an SVD.
}
  \item{left.additions}{
  A list of additions to the ggplots of the left singular vectors.  Can be used to set aesthetics such as color, etc.
}
  \item{right.additions}{
  A list of additions to the ggplots of the right singular vectors.  Can be used to set aesthetics such as color, etc.
}
  \item{factor.labels}{
  The factor labels.
}
}

\details{
  Plots of the left singular vectors are shown on the left, and plots of the right singular vectors are shown on the right.  The diagonal shows squares with side lengths proportional to the singular values.  
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


