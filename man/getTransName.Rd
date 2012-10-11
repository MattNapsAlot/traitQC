\name{getTransName}
\alias{getTransName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Transformation Name }
\description{
  Get the transformation name corresponding to the input lambda value.
}
\usage{
getTransName(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ lambda value }
}

\value{
  Text representation of the transformation specified by the value of lambda.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{getLambda}}, \code{\link{boxcox}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x){
	lambdas <- c(-2, -1, -0.5, 0 , 0.5, 1, 2);
	nms <- c('inverse_square','inverse','inverse_sqrt','log','sqrt','normal','square');
	if(!(x \%in\% lambdas)) stop(paste('invalid transform:',x));
	ind <- which(x == lambdas);
	return(nms[ind]);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lambda }
\keyword{ box }% __ONLY ONE__ keyword per line
\keyword{ cox }
\keyword{ power }
\keyword{ transform }