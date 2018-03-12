\name{collapse.replicates}
\alias{collapse.replicates}

\title{
Collapse Replicates
}

\description{
This function is intended for use in conjunction with \code{\link{RUVIII}}, specifically when using the \code{average=TRUE} option.  When using the \code{average=TRUE} option, the adjusted data matrix has only one row for each set of replicates.  In other words, each set of replicate rows in the orginal data matrix is "collapsed" into a single row in the adjusted data matrix.  This function similarly collapses the rows of a dataframe of covariates.  Only covariates that are constant within each set of replicates are retained.    
}

\usage{
collapse.replicates(df, M)  
}

\arguments{
  \item{df}{
   A dataframe.
}
  \item{M}{
   The replicate structure.  See \code{\link{RUVIII}} for details.
}

}

\value{
A sub-dataframe of df.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}


\seealso{
\code{\link{RUVIII}}
}
