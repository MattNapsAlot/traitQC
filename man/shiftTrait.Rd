\name{shiftTrait}
\alias{shiftTrait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Shift Trait Values }
\description{
  Add a constant value to the input values so that the input array contains no zero or negative values.
}
\usage{
shiftTrait(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input values }
}
\details{
  If x contains any values <= 0, compute a shift constant 
  and adds it to the input values.
  \deqn{shift=abs(min(x))+1}
  \deqn{y=x+shift}
  
  If all values in x are > 0, then \eqn{shift=0} and \eqn{y=x}
}
\value{
  A list with the following components:
  \item{shift }{Shift applied to the input values}
  \item{ y }{The shifted values}
}
\author{ Matt Furia }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(y){
	ret = list(y=y,shift=0);
	if(any(y <= 0)){
		ret$shift = abs(min(y)) + 1;
		ret$y = y + ret$shift;
	}
	return(ret);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ shift }
