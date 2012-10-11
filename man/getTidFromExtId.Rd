\name{getTidFromExtId}
\alias{getTidFromExtId}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get TIDs from DB }
\description{
  Query MegaX_Works to retrieve TIDs.
}
\usage{
getTidFromExtId(ids, na.rm = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{ External names }
  \item{na.rm}{ boolean indicating whether individuals missing from the DB should be ommitted from the return value }
}
\value{
  A list containing the external_name/tid map:
  \item{external_name }{external name}
  \item{tid }{tid}
  ...
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Matt Furia }
\note{ ~~further notes~~ 
  Currently this function is only supported for the Windows OS due to it's
  dependence on the RODBC package.  Additionally, an ODBC data sorce must
  be set up which points to the server where the MegaX_works DB resides.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(ids,na.rm=F){
	# get the tid/external ID map from db
	require(RODBC);
	conn = odbcConnect("USSEDB1604");
	stmt = "select distinct tid, external_name, cross_name from individual where external_name in(";
	stmt = paste(stmt,paste(paste("'",ids,"'",sep=""),collapse=','),")",sep="");

	map = sqlQuery(conn, stmt);
	ind = match(ids,map$external_name);
	odbcClose(conn);
	if(na.rm){
		ind = ind[!is.na(ind)];
	}
	return(map$tid[ind]);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tid }
