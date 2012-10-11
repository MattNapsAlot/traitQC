\name{pivotTable}
\alias{pivotTable}
\alias{unpivotTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Table Pivoting }
\description{
  Convert data frames between tabular and attribute/value formats.
}
\usage{
pivotTable(t, pivot.fields, na.rm=F)
unpivotTable(t, index.fields=c('id','interval'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Data frame }
  \item{pivot.fields}{ Fields to be used for pivoting }
  \item{index.fields}{ Index fields to be used for un-pivoting }
  \item{na.rm}{boolean indicating whether NA attribute values should be removed from pivoted table}
}
\details{
  When pivoting, the new data frame frame will have the fields 'trait_name' and 'trait_value' in addition to the
  \code{pivot.fields}.  
  
  When un-pivoting, the input data frame must contain the fields 'trait_name' and 'trait_value' and
  the new data frame will contain the \code{index.fields}plus a field for each unique value 
  contained in the 'trait_name' field of the input data frame.
}
\value{
  Data frame containing data in the new format. 
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }

\seealso{ \code{\link{unpivotTable}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pivot }
