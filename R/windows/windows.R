#
#get the tid from the external id;
#
getTidFromExtId = function(ids,na.rm=F){
	# get the tid/external ID map from db
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