\name{ruv_cancorplot}
\alias{ruv_cancorplot}

\title{
RUV Canonical Correlation Plot
}

\description{
Canonical correlation plot
}

\usage{
ruv_cancorplot(Y, X, ctl, W1 = NULL, W2 = NULL)
}

\arguments{
  \item{Y}{
  The data matrix.  Rows are observations and columns are features (e.g. genes).  
}
  \item{X}{
  Factor(s) of interest.  Can be a vector, factor, matrix, or dataframe.  Must have the same length (or number of rows) as the number of row of Y.  
}
  \item{ctl}{
  Index of negative controls.
}
  \item{W1}{
  Optional.  The left singular vectors of \code{Y}.  Can be included to speed up execution.
}
  \item{W2}{
  Optional.  The left singular vectors of \code{Y[,ctl]}.  Can be included to speed up execution.
}
}

\details{
  Plots, as a function of k, the square of of the first canonical correlation of \code{X} and the first k left singular vectors of \code{Y} (and also, similarly, \code{Y[,ctl]}).
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


