\name{getLambda}
\alias{getLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Lambda Value }
\description{
  Return the value of lambda that is closest to the input value.
}
\usage{
getLambda(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{ lambda estimate from boxcox }
}
\details{
  Possible values of lambda are {-2, -1, -0.5, 0, 0.5, 1, 2}.  This function
  returns the value from this set that is the closest to the input value.
}
\value{
  lambda.
}
\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{boxcox}}, \code{\link{bxcxTrans}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(input){
	lambdas <- c(-2, -1, -0.5, 0 , 0.5, 1, 2);
	lambda <- c();
	for(x in input){
		dif <- x - lambdas;
		ind <- which.min(dif^2);
		if(length(ind) > 1){
			if(sign(x) < 0){
				l <- lambdas[ind[1]];
			}else{
				l <- lambdas[ind[2]];
			}
		}else{
			l <- lambdas[ind];
		}
		lambda <- c(lambda,l);
	}
	return(lambda);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lambda }
\keyword{ box }
\keyword{ cox }
