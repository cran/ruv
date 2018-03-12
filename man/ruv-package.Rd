\name{ruv-package}
\alias{ruv-package}
\alias{ruv}
\docType{package}

\title{
Detect and Remove Unwanted Variation using Negative Controls 
}

\description{
Implements the 'RUV' (Remove Unwanted Variation) algorithms.  These algorithms attempt to adjust for systematic errors of unknown origin in high-dimensional data.  The algorithms were originally developed for use with genomic data, especially microarray data, but may be useful with other types of high-dimensional data as well.  These algorithms were proposed by Gagnon-Bartsch and Speed (2012), and by Gagnon-Bartsch, Jacob and Speed (2013).  The algorithms require the user to specify a set of negative control variables, as described in the references.  The algorithms included in this package are 'RUV-2', 'RUV-4', 'RUV-inv', 'RUV-rinv', 'RUV-I', and RUV-III', along with various supporting algorithms.  
}

\details{
\tabular{ll}{
Package: \tab ruv\cr
Type: \tab Package\cr
Version: \tab 0.9.7\cr
Date: \tab 2018-03-12\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
URL: \tab \href{http://www-personal.umich.edu/~johanngb/ruv/}{http://www-personal.umich.edu/~johanngb/ruv/} \cr
}
}

\keyword{models} 
\keyword{multivariate} 
 
\seealso{ 
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{variance_adjust}}, 
\code{\link{RUVI}}, \code{\link{RUVIII}} 
} 
 


\author{ 
  Johann Gagnon-Bartsch <johanngb@umich.edu> 
} 
 
\references{ 
Using control genes to correct for unwanted variation in microarray data.  Gagnon-Bartsch and Speed, 2012.  Available at: \href{http://biostatistics.oxfordjournals.org/content/13/3/539.full}{http://biostatistics.oxfordjournals.org/content/13/3/539.full}.   
 
Removing Unwanted Variation from High Dimensional Data with Negative Controls.  Gagnon-Bartsch, Jacob, and Speed, 2013.  Available at: \href{http://statistics.berkeley.edu/tech-reports/820}{http://statistics.berkeley.edu/tech-reports/820}.  
} 
 
\note{ 
Additional resources can be found at \href{http://www-personal.umich.edu/~johanngb/ruv/}{http://www-personal.umich.edu/~johanngb/ruv/}.   
} 
