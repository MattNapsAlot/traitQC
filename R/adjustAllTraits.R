# TODO: Add comment
# 
# Author: furia
###############################################################################


adjustAllTraits = function(trait,covar,include.unadjusted=F,relevel.covar=T,thr=0.01,marg.thr=0.01,id.fname='external_name'){
	require(traitQC)
	if(any(is.na(covar))) stop('ERROR: Missing some covariate values')
	if(relevel.covar){
		#relevel covariates
		mk = covar[,id.fname] %in% unique(trait[,id.fname])
		covar = covar[mk,]
		for(n in setdiff(names(covar),id.fname)){
			covar[,n] = covarAdj(covar[,n]) 
		}
	}
	trait.final = trait[0,]
	for(tr in unique(trait$trait_name)){
		mk.tr = trait$trait_name == tr
		for(tp in unique(trait$timepoint[mk])){
			mk = mk.tr & trait$timepoint == tp
			tmp = trait[mk,]
			trait.final = rbind(trait.final,adjustTrait(tmp,covar,thr=thr,marg.thr=marg.thr,id.fname=id.fname))
			rm(tmp)
		}
	}
	if(include.unadjusted){
		trait.final = rbind(trait,trait.final)
	}
	return(trait.final)
}

applyAdjust = function(trait,covar,model){
	d = merge(trait,covar)
	fit = lm(model,data=d)
	return(array(fit$resid) + fit$coefficients['(Intercept)'])
}

dropCovars = function(trait,covar,drop.covar){
	for(tr in unique(trait$trait_name)){
		mk.tr = trait$trait_name == tr
		for(tp in unique(trait$timepoint[mk])){
			mk = mk.tr & trait$timepoint == tp
			tmp = trait[mk,]
			trait.final = rbind(trait.final,adjustTrait(tmp,covar,thr=thr,marg.thr=marg.thr,id.fname=id.fname))
			rm(tmp)
		}
	}
}


