\name{ruv_hist}
\alias{ruv_hist}

\title{
RUV P-value Histogram Plot
}

\description{
Plots a histogram of p-values returned by a call to \code{\link{ruv_summary}}
}

\usage{
ruv_hist(fit, X.col = "all", breaks = c(0, 0.001, 0.01, 0.05, seq(0.1, 1, by = 0.1)))
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{X.col}{
  Which column of the X matrix to make the plot for, i.e. which factor's p-values to plot.  Can be either an integer or a character string.  Or, if "all" (the default), use the F-test p-values.
}
  \item{breaks}{
  Breakpoints of the histogram.
}
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


