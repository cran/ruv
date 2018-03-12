\name{design.matrix}
\alias{design.matrix}

\title{
Design Matrix
}

\description{
Creates a design matrix.
}

\usage{
design.matrix(a, name = "X", remove.collinear = TRUE, include.intercept = TRUE)
}

\arguments{
  \item{a}{
    Object from which to create a design matrix.  Can be a vector, matrix, factor, or dataframe.
}
  \item{name}{
    Name of the design matrix.  Used to name columns that aren't already named (e.g. X1, X2, etc.)
}
  \item{remove.collinear}{
    Will remove columns that are collinear, to ensure the design matrix is full rank.
}
  \item{include.intercept}{
    Add an intercept to the matrix if one is not included already.
}
}

\details{
    Numerical vectors are not modified.  Factors are converted to dummy variables.  Character vectors are converted to factors, and then to dummy variables.
}

\value{
A matrix.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}
