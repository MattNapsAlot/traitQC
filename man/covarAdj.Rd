\name{covarAdj}
\alias{covarAdj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adjust/Re-level Covariates }
\description{
  Adjust or re-level the covariate depending on whether it is a factor or a continuous variable.
}
\usage{
covarAdj(covar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covar}{ An array of covariate values. Values can be either factors or continuous variables.  }
}
\details{
  For continuous covariates, the covariate values are adjusted by subtracting
  the mean from each value.  If the covariate is a factor, it is re-leveled using
  the \code{relevel} function from the \code{stats} package. The reference level 
  for re-leveling is the one that occurrs most frequently in the input data.  
}
\value{
  An array of adjusted or re-leveled covariate values.
}
#\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{relevel}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(covar){
	if(is.character(covar)) covar = factor(covar);
	if(is.numeric(covar)){
		covar = covar - mean(covar, na.rm =T);
	}else if(is.factor(covar)){
		ref = names(summary(covar))[which(summary(covar)==max(summary(covar)))][1];
		covar = relevel(covar,ref);
	}
	return(covar);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ covariate }
\keyword{ adjust }
\keyword{ relevel }% __ONLY ONE__ keyword per line
