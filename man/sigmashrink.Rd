\name{sigmashrink}
\alias{sigmashrink}

\title{
Empirical Bayes shrinkage estimate of sigma^2
}

\description{
This function (re)implements the empirical bayes shrinkage estimate of Smyth (2004), which is also implemented in the Limma package.  This function is normally called from the function \code{\link{variance_adjust}}, and is not normally intended for stand-alone use.
}

\usage{
sigmashrink(s2, d)
}

\arguments{
  \item{s2}{
  "Standard" estimates of sigma^2
}
  \item{d}{
  "Standard" degrees of freedom of the residuals
}
}

\value{
A list containing
  \item{sigma2}{Estimates of sigma^2 using the empirical bayes shrinkage method of Smyth (2004)}
  \item{df}{Estimate of degrees of freedom using the empirical bayes shrinkage method of Smyth (2004)}
}


\references{
Linear models and empirical bayes methods for assessing differential expression in microarray experiments.  Smyth, 2004.

Using control genes to correct for unwanted variation in microarray data.  Gagnon-Bartsch and Speed, 2012.  Available at: \href{http://biostatistics.oxfordjournals.org/content/13/3/539.full}{http://biostatistics.oxfordjournals.org/content/13/3/539.full}. 

Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}

\seealso{
\code{\link{variance_adjust}}
}

