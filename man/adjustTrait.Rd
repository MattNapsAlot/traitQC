\name{adjustTrait}
\alias{adjustTrait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adjust and Transform Clinical Traits For Covariates}
\description{
  Performs standard QC on a clinical trait by adjusting for significant covariates, transforming for normality and shifting values to eliminate non-positive values.
}
\usage{
adjustTrait(t, covar, thr = 0.01, marg.thr = 0.01, id.fname = "external_name")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Data frame containing the trait data. See Details. }
  \item{covar}{ Covariates data frame. See Details.}
  \item{thr}{ significance threshold for the final adjustment model. }
  \item{marg.thr}{ significance threshold for assesing marginal effects. }
  \item{id.fname}{ fieldname in \code{t} and \code{covar} that contains the individual identifier }
}
\details{
  In addition to a field containing the individual identifier, which must be named \code{id.fname}, the data frame 
  containing the trait data must have the following fields:
	trait_name: The trait name for the observation
	trait_value: The value

  All fields in the covariates data frame other than the field named \code{id.fname} will be considered covariates
  and will considered for inclusion in the adjustment model.
  
  Prior to model selection, the covariates are adjusted using \code{covarAdj} and
  the trait is shifted using \code{shiftTrait}.

  Model selection begins by checking the marginal significance of each
  covariate against the trait.  A multivariate model is then build using the
  covariates that meet the significance threshold \code{marg.thr}.  \code{stepAIC}
  is used to determine the final model and covariates that do not meet the 
  significance threshold \code{thr} are removed from the model.
  
  \code{boxcox} is used to determine the appropriate power transform.  If necessary,
  the trait is transformed, the model re-fit and effects that don't meet the significance 
  threshold \code{thr} are removed.
  
  Adjusted trait values are computed as the residuals from the model plus the intercept.
}

\value{
  Data frame containing the new traits.  The returned trait values will potentially be shifted to eliminate 
  non-positive values, adjusted for covariates and transformed for normality.  The data frame will contain
  intermediate values in addition to the fully transformed and adjusted values.  For example, if it is determined
  that a trait needs to be adjusted for covariates and transformed, the data frame will include the adjusted but
  untransformed values and the transformed but unadjusted values in addition to the fully transformed/adjusted values.

  In addition to the field containing the individual identifiers (named \code{id.fname}), the data frame will contain
  the following fields:
  timepoint: The timepoint for the observation
  trait_name: The trait name
  shift: The constant value added to the original trait value prior to performing any adjustments/transformations
  adjustment: The covariates used in the adjustment model.
  trait_value: The new trait value after applying the shift, transformation and adjustment.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }
%\note{ ~~further notes~~ 
%~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{\link{stepAIC}}, \code{\link{covarAdj}}, \code{\link{shiftTrait}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ clinical }
\keyword{ trait }% __ONLY ONE__ keyword per line
\keyword{ QC }

