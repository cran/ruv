\name{ruv_varianceplot}
\alias{ruv_varianceplot}

\title{
RUV Variance Plot
}

\description{
  A scatter plot of (squared) coefficient estimates against variance estimates.
}

\usage{
ruv_varianceplot(fit, X.col = 1, power = 1/4)
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{X.col}{
  Which column of the X matrix to make the plot for.  Can be either an integer or a character string.
}
  \item{power}{
  Power transformation of the x and y axes.  Default is fourth root.
}
}

\details{
  A black curve is also plotted, showing the estimated variances of the coefficient estimates.
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


