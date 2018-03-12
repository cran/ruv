\name{ruv_summary}
\alias{ruv_summary}

\title{
RUV Summary
}

\description{
Post-process and summarize the results of call to RUV2, RUV4, RUVinv, or RUVrinv.
}

\usage{
ruv_summary(Y, fit, rowinfo=NULL, colinfo=NULL, colsubset=NULL, sort.by="F.p", 
            var.type=c("ebayes", "standard", "pooled"),
            p.type=c("standard", "rsvar", "evar"), min.p.cutoff=10e-25)
}

\arguments{
  \item{Y}{
  The original data matrix used in the call to RUV2/4/inv/rinv
}
  \item{fit}{
  A RUV model fit (a list), as returned by RUV2 / RUV4 / RUVinv / RUVrinv
}
  \item{rowinfo}{
  A matrix or dataframe containing information about the rows (samples).  This information is included in the summary that is returned.
}
  \item{colinfo}{
  A matrix or dataframe containing information about the columns (features, e.g. genes).  This information is included in the summary that is returned.
}
  \item{colsubset}{
  A vector indexing the features of interest.  Only only data on these features will be returned.
}
  \item{sort.by}{
  An index variable; which column of \code{C} (see below) should be used to sort the features.  The default is "F.p", meaning that features will be sorted by the F-test p-value.  If \code{NULL}, features will note be sorted.
}
  \item{var.type}{
  Which type of estimate for sigma2 should be used from the call to \link{variance_adjust}?  The options are "ebayes", "standard", or "pooled."  See \link{variance_adjust} for details.
}
  \item{p.type}{
  Which type of p-values should be used from the call to \link{variance_adjust}?  The options are "standard", "rsvar", or "evar". 
}
  \item{min.p.cutoff}{
  p-values below this value will be changed and set equal to this value.  Useful for plotting p-values on a log scale.
}
}

\details{
This function post-processes the results of a call to RUV2/4/inv/rinv and then nicely summarizes the output.  The post-processing step primarily consists of a call to \link{variance_adjust}, which computes various adjustments to variances, t-statistics, and and p-values.  See \link{variance_adjust} for details.  The \code{var.type} and \code{p.type} options determine which of these adjustments are used.  An additional post-processing step is that the column means of the \code{Y} matrix are computed, both before and after the call to \code{RUV1} (if \code{eta} was specified).  

After post-processing, the results are summarized into a list containing 4 objects: 1) the data matrix \code{Y}; 2) a dataframe \code{R} containing information about the rows (samples); 3) a dataframe \code{C} containing information about the columns (features, e.g. genes), and 4) a list \code{misc} of other information returned by RUV2/4/inv/rinv. 

Finally, if \code{colsubset} is specified, then \code{C} is subset to include only the features of interest (as are the relevant entries of \code{misc} that are used to compute projection plots).  If \code{sort.by} is specified, the features will also be sorted.
}

\value{
A list containing:
    \item{Y}{The original data matrix.}
    \item{R}{A dataframe of row-wise information, including \code{X}, \code{Z}, and any other data passed in with \code{rowinfo}}
    \item{C}{A dataframe of column-wise information, including p-values, estimated regression coefficients, estimated variances, column means, an index of the negative controls, and any other data passed in with \code{colinfo}.}
    \item{misc}{A list of additional information returned by RUV2/4/inv/rinv}
}

\author{
Johann Gagnon-Bartsch
}

\seealso{
\code{\link{RUV2}}, \code{\link{RUV4}}, \code{\link{RUVinv}}, \code{\link{RUVrinv}}, \code{\link{variance_adjust}}
}

