AddPvalueToHDF5<-function(edc,files,idType="tid",hc,gopts,h5.paths,verbose=T)
{
	block.size = 1000
	default.paths = list(
			indiv.path = '/all.data/individuals',
			exprid.path = '/all.data/repIDs',
			pval.path = '/all.data/pvalue',
			expr.path = '/all.data/repData'
	)
	require(bat)
	memory.size(3996)
	
	#set up the h5 path defaults
	if(missing(h5.paths)){
		h5.paths = default.paths
	}
	for(p in setdiff(names(default.paths),names(h5.paths))){
		h5.paths[p] = default.paths[p];
	}
	
	#get the pvalues from resolver
	if(missing(hc)){
		hc = ResolverLogin("infmse", "bat")
	}
	if(missing(gopts)){
		gopts<-list(OrderBy="none",IncludeFlaggedData="False")
	}
	if(verbose) print(paste("Getting EG for EDC(s): ", edc))
	eg<-GetResolverEG(hc,c("pvalue"),"mouse",RDAids("edc",edc),lod="reporter",tod="re",bsa=TRUE,bsa.simplify=FALSE,options=gopts,debug=FALSE)
	
	#rbinary file name. for debugging use
	#fname = paste(gsub(';','_',edc),'rbin',sep='.');
	
	#HACK this needs to be fixed and will cause incorrect pvalues. Need a way to match hdf files to chip barcodes
	#drop duplicates
	mk = table(eg$bsainfo.exp[[idType]]) > 1
	if(sum(mk > 0)){
		dup.ids = names(table(eg$bsainfo.exp[[idType]]))[mk]
		for(id in dup.ids){
			ind = which(eg$bsainfo.exp[[idType]] == id)
			#keep the first one
			ind = ind[-1]
			eg$ced$pvalue = eg$ced$pvalue[-ind,]
			eg$bsainfo.exp[[idType]] = eg$bsainfo.exp[[idType]][-ind]
			print(paste("WARNING: found duplicate values for tid", id, "pvalues for this tid should be treated with skepticism"));
		}
	}
	#END HACK

	for (f in files)
	{
		print(paste("Adding pvalues to ", f))
		x <- NULL
		x$individuals <- HDF5ReadData(f, h5.paths$indiv.path)
		x$repIDs  <- HDF5ReadData(f, h5.paths$exprid.path)
		#pre make the size of the data
		x$pvalue = matrix(nrow=length(x$repIDs),ncol=length(x$individuals))
		x$pvalue[]<-NA
		
		startindex = 0
		write.cnt = 1
		start = startindex + 1
		end = min(write.cnt*block.size, dim(x$pvalue)[1])
       
		
		#match the individuals
		if(idType=="tid")
		{
			ids = intersect(eg$bsainfo.exp[[idType]],x$individuals)
			mk.eg = eg$bsainfo.exp[[idType]] %in% ids
			mk.x = x$individuals %in% ids
			ind<-match(x$individuals[mk.x],eg$bsainfo.exp[[idType]][mk.eg])
		}
		else
		{
			stop("unsupported id type")
			ind<-match(eg$bsainfo.exp[[idType]]$data,x$individuals)
		}
		mk1<-eg$lodids %in% x$repIDs
		mk2<-x$repIDs %in% eg$lodids
		ind1<-match(x$repIDs[mk2],eg$lodids[mk1])

		#place pvalue into the rbin
		x$pvalue[mk2,mk.x]<- t(eg$ced$pvalue[mk.eg,mk1][ind,ind1])
		while(end != dim(x$pvalue)[1]){
			end = min(write.cnt*block.size, dim(x$pvalue)[1])
			#save the pvalue in the hdf5 value
			if(startindex == 0){
				HDF5WriteData(f,h5.paths$pval.path,x$pvalue[start:end,],list(overwrite=T))
			}else{
	        	HDF5WriteData(f,h5.paths$pval.path,x$pvalue[start:end,],list(overwrite=T,startindex=startindex))
			}
	        startindex = end
			start = startindex + 1
			write.cnt = write.cnt + 1
		}
		rm(x)
	}
	rm(eg)
}


