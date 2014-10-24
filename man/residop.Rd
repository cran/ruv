\name{residop}
\alias{residop}

\title{
Residual Operator
}

\description{
Applies the residual operator of a matrix B to a matrix A.
}

\usage{
residop(A, B)
}

\arguments{
  \item{A}{
     A matrix with m rows.
}
  \item{B}{
     Another matrix with m rows.
}
}

\details{
The columns of B must be linearly independent.
}

\value{
The matrix A after projection into the orthogonal complement of the column space of B.
}




