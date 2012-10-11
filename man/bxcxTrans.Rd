\name{bxcxTrans}
\alias{bxcxTrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Power Transform}
\description{
  Apply the power transform to a response variable.
}
\usage{
bxcxTrans(y, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ The values to be transformed }
  \item{lambda}{ \code{lambda} value for the transform. }
}
\details{
  Apply the power transform to the input values and return the result. Valid values
  of \code{lambda} are {-2, -1, -0.5, 0, 0.5, 1, 2}.  for \code{lambda}=0 the input
  values are transformed using log(\code{y}), otherwise the values are transformed
  using \eqn{y^lambda}. Input values must be greater than 0.
}
\value{
  An array of transformed values.
}
\references{ \url{http://en.wikipedia.org/wiki/Power_transform} }
\author{ Matt Furia }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{getLambda}}, \code{\link{getTransName}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(y, lambda){
	lambdas = c(-2, -1, -0.5, 0 , 0.5, 1, 2);
	if(length(lambda) != 1) stop("invalid input for lambda");
	if(!(lambda \%in\% lambdas)) stop(c("invalid lambda value:", lambda));
	if(lambda == 0){
		ret = do.call("log",list(y));
	}else{
		ret = do.call("power",list(y,lambda));
	}
	return(ret);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ box }
\keyword{ cox }% __ONLY ONE__ keyword per line
\keyword{ transform }
\keyword{ power }
