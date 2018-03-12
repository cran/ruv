\name{replicate.matrix}
\alias{replicate.matrix}

\title{
Replicate (Mapping) Matrix
}

\description{
For use with \code{\link{RUVIII}}, generates a mapping matrix that describes the replicate structure.  
}

\usage{
replicate.matrix(a, burst=NULL, return.factor=FALSE, name="M", sep="_", burstsep = "_")
}

\arguments{
  \item{a}{
   An object that describes the replicate structure.  Can be a vector, matrix, factor, or dataframe.  If a vector, it is converted to a factor.  If a factor, each level of the factor is taken to represent a set of replicates.  If a matrix: First it is determined whether it is already a mapping matrix; if so, the matrix is returned unchanged; if not, the matrix is converted to a dataframe.   If a dataframe: Each column is converted to a factor.  A new factor is then created with levels for every possible combination of factor levels in the dataframe.  For example, if the dataframe contains two factors, patientID and sampleDate, the new factor will have a unique level for each (observed) combination of patientID and sampleDate.  Thus observations will be considered replicates if they have identical values for BOTH patientID and sampleDate.
}
  \item{burst}{
   A character vector, containing the names of factor levels to be "burst."  When a factor level is burst, the corresponding observations are no longer replicates; they become singletons.
}
  \item{return.factor}{
   Return a factor instead of the mapping matrix.  This may be useful in two situations: (1) When the input is a mapping matrix, and it is desired to convert it back to a factor; (2) When making repeated calls to \code{replicate.matrix} in order to define the replicates in several steps.  Example of (2): Suppose there are 4 patients and 3 sample dates.  We wish to designate as replicates observations that have the same patient ID and sample date, but only for the first two sample dates; none of the observations in the third sample date should be considered replicates.  We would first call \code{replicate.matrix} using only the sampleDate factor, bursting the third level, and returning another factor.  We would then call \code{replicate.matrix} again, this time with a dataframe containing patientID and the bursted sampleDate.  See below for example code.
}
  \item{name}{
   Name of the mapping matrix.  Used to name columns that aren't already named (e.g. M1, M2, etc.)
}
  \item{sep}{
   Text separating the level names of different factors when they are combined.  
}
  \item{burstsep}{
   Text appended to factor level names when bursting a factor.  This text is then followed by a number.  Example: if the factor level to be burst is "June29", and burstsep is the default value of "_", then the new levels will be "June29_1", "June29_2", etc.
}
}

\section{Warning }{Be sure to change the default values of \code{sep} and \code{burstsep} if there is any risk of factor level naming conflicts (e.g. if existing factors already have level names like "patient_1", "patient_2", etc.}

\value{
A matrix or a factor, depending on the value of \code{return.factor}.
}

\author{
Johann Gagnon-Bartsch \email{johanngb@umich.edu}
}

\examples{
# Define patientID and sampleDate
patientID = paste("patient", rep(1:4, each=6), sep="")
#print(patientID)
sampleDate = paste("June", rep(c(12,17,29), 8), sep="")
#print(sampleDate)
# Create a mapping matrix, where every unique 
# patientID / sampleDate combination define a set of replicates
M = replicate.matrix(data.frame(patientID, sampleDate))
#print(M)
# Convert M back to a factor
M = replicate.matrix(M, return.factor=TRUE)
#print(M)
# Create a factor for sampleDate, but burst the third date
temp = replicate.matrix(sampleDate, burst="June29", return.factor=TRUE)
#print(temp)
# Create a mapping matrix as described above in the description of return.factor
M = replicate.matrix(data.frame(temp, patientID))
#print(M)
}

\seealso{
\code{\link{RUVIII}}
}
