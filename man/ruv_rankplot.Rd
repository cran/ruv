\name{ruv_rankplot}
\alias{ruv_rankplot}

\title{
RUV Rank Plot
}

\description{
A plot showing the number of positive controls to be found within the N top-ranked features, as a function of N.  The ranking of the features is by p-value. 
}

\usage{
ruv_rankplot(fit, pctl, X.col = "all", uniform.lines = 0)
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{pctl}{
  Either an integer or character string specifying which column of \code{fit$C} to be used as positive controls.  (Must be a logical vector).  Alternatively, may some other index vector specifying the positive controls; importantly, in this case, the index vector must index the features as they are sorted in \code{fit$C}.
}
  \item{X.col}{
  Which column of the X matrix to make the plot for.  Can be either an integer or a character string.  Or, if "all" (the default), use the F-test p-values.
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


