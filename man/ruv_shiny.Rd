\name{ruv_shiny}
\alias{ruv_shiny}

\title{
RUV Shiny App
}

\description{
  A Shiny App that allows quick exploration of a dataset using RUV methods.
}

\usage{
ruv_shiny(Y, rowinfo, colinfo, options = list(port = 3840)) 
}

\arguments{
  \item{Y}{
  The data matrix.  Rows are observations and columns are features (e.g. genes).
}
  \item{rowinfo}{
  A dataframe of information about the observations.  Should have the same number of rows as Y.  Should contain at least one column that can be used as either a factor of interest or to define replicates.
}
  \item{colinfo}{
  A dataframe of information about the observations.  Should have the same number of rows as Y.  Should contain at least one column that is a logical vector that can be used to define negative controls.
}
  \item{options}{
  A list of options to pass to the shinyApp function.
}
}

\value{
None.  Calls shinyApp.
}

\author{
Johann Gagnon-Bartsch
}


