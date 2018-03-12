\name{ruv_projectionplot}
\alias{ruv_projectionplot}

\title{
RUV Projection Plot
}

\description{
Projection plot of an RUV regression fit.
}

\usage{
ruv_projectionplot(fit, X.col = 1, factor = "gradient", adjusted = TRUE)
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{X.col}{
  Which column of the X matrix to make the plot for.  Can be either an integer or a character string.
}
  \item{factor}{
  Which unwanted factor to use (horizontal axis).  Must be either an integer or the character string "gradient".
}
  \item{adjusted}{
  Whether the plot should be adjusted for unwanted factors other than the one being plotted.  Not relevant when \code{factor = "gradient"}.
}
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


