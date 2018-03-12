\name{ruv_ecdf}
\alias{ruv_ecdf}

\title{
RUV P-value Empirical CDF Plot
}

\description{
Plots an ECDF of p-values returned by a call to \code{\link{ruv_summary}}
}

\usage{
ruv_ecdf(fit, X.col = "all", power = 1, uniform.lines = 0)
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{X.col}{
  Which column of the X matrix to make the plot for, i.e. which factor's p-values to plot.  Can be either an integer or a character string.  Or, if "all" (the default), use the F-test p-values.
}
  \item{power}{
  A power transformation of the x and y axes.  For example, set to 1/2 for a square-root transformation.  This can help to see the behavior of the ECDF near 0.
}
  \item{uniform.lines}{
  A vector of values between 0 and 1, or NULL.  If specified, light gray lines will be drawn, showing (locally) uniform distributions. 
}
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


