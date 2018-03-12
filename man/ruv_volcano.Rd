\name{ruv_volcano}
\alias{ruv_volcano}

\title{
RUV Volcano Plot
}

\description{
  A scatter plot of negative log p-values against coefficient estimates, commonly known as a volcano plot
}

\usage{
ruv_volcano(fit, X.col = 1)
}

\arguments{
  \item{fit}{
  The results of a call to \code{\link{ruv_summary}}.
}
  \item{X.col}{
  Which column of the X matrix to make the plot for.  Can be either an integer or a character string.
}
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}


