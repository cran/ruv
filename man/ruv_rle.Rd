\name{ruv_rle}
\alias{ruv_rle}

\title{
RUV RLE Plot
}

\description{
An RLE (Relative Log Expression) Plot
}

\usage{
ruv_rle(Y, rowinfo = NULL, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), ylim = c(-0.5, 0.5))
}

\arguments{
  \item{Y}{
  The data matrix.  Rows are observations and columns are features (e.g. genes).
}
  \item{rowinfo}{
  A dataframe of information about the observations.  Should have the same number of rows as Y.  This information will be included in the ggplot, and can be used for setting aesthetics such as color.
}
  \item{probs}{
  The percentiles used to construct the boxplots.  By default, whiskers are drawn to the 5th and 95th percentiles.  Note that this is non-standard for boxplots.
}
  \item{ylim}{
  Limits of the y axis.  Defaults to (-0.5, 0.5) so that the plots are always on the same scale and can be easily compared.
}
}

\value{
A ggplot.
}

\author{
Johann Gagnon-Bartsch
}

\references{
Irizarry, R. A.,  Bolstad, B. M., Collin, F., Cope, L. M., Hobbs, B., and Speed, T. P. (2003).  Summaries of Affymetrix GeneChip probe level data.  Nucleic acids research, 31(4):e15.
}
