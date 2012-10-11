\name{pwr}
\alias{pwr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exponent }
\description{
  Raise the input values to the power specified
}
\usage{
pwr(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Starting values }
  \item{y}{ The exponent(s) }
}
\details{
  If length(y) == 1, all the x's will be raised the the same power.  If 
  length(y) == length(x), each x will be raised the the corresponding power
  in y.
}
\value{
  The computed exponentials.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ exponent }
\keyword{ power }% __ONLY ONE__ keyword per line
