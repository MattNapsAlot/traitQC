
getLambda <- function(input){
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

getTransName <- function(x){
	lambdas <- c(-2, -1, -0.5, 0 , 0.5, 1, 2);
	nms <- c('inverse_square','inverse','inverse_sqrt','log','sqrt','normal','square');
	if(!(x %in% lambdas)) stop(paste('invalid transform:',x));
	ind <- which(x == lambdas);
	return(nms[ind]);
}



pwr = function(x,y){
	if(length(y) != 1 & length(y) != length(x)) stop('unexpected input for y');
	return (x^y);
}


bxcxTrans = function(y, lambda){
	lambdas = c(-2, -1, -0.5, 0 , 0.5, 1, 2);
	if(length(lambda) != 1) stop("invalid input for lambda");
	if(!(lambda %in% lambdas)) stop(c("invalid lambda value:", lambda));
	if(any(y <= 0)) stop("input values must be greater than 0.");
	if(lambda == 0){
		ret = do.call("log",list(y));
	}else{
		ret = do.call("pwr",list(y,lambda));
	}
	return(ret);
}

covarAdj = function(covar){
	if(is.character(covar)) covar = factor(covar);
	if(is.numeric(covar)){
		covar = covar - mean(covar, na.rm =T);
	}else if(is.factor(covar)){
		ref = names(summary(covar))[which(summary(covar)==max(summary(covar)))][1];
		covar = relevel(covar,ref);
	}
	return(covar);
}

pivotTable = function(t,pivot.fields,na.rm=F){
	if(!is.data.frame(t)){
		stop("input must be a data frame");
	}
	
	pivot.fields = unique(pivot.fields);
	if(!all(pivot.fields %in% names(t))){
		stop(paste("invalid pivot fields:",pivot.fields));
	}
	
	tmp = t[pivot.fields];
	tmp$trait_name = NA;
	tmp$trait_value = NA;
	ret = tmp[0,];
	nms = names(t)[!(names(t) %in% pivot.fields)];
	
	for(fn in nms){
		tmp$trait_name = fn;
		tmp$trait_value = t[,fn];
		ret = rbind(ret,tmp);
	}
	if(na.rm){
		mk = is.na(ret$trait_value);
		ret = ret[!mk,];
	}
	return(ret);
}

unpivotTable = function(t,index.fields=c('id','interval')){
	if(!all(index.fields %in% names(t))) stop("invalid index fields");
	fnames = c('trait_value','trait_name');
	if(!all(fnames %in% names(t))) stop("Input data frame must contain the fields 'trait_name' and 'trait_value'");
	ret = unique(t[index.fields]);
	rownames(ret) = 1:dim(ret)[1];
	tnames = unique(as.character(t$trait_name));
		
	for(tr in tnames){
		mk = as.character(t$trait_name) == tr;
		tmp = unique(t[mk,c(index.fields,'trait_value')]);
		names(tmp)[names(tmp) == 'trait_value'] = tr; 
		if(dim(tmp)[1] > dim(ret)[1]) stop("found more than one value for a single trait for an individual");
		ret = merge(ret,tmp,all.x=T);
	}
	return(ret);
}


shiftTrait = function(x){
	ret = list(y=x,shift=0);
	if(any(x <= 0)){
		ret$shift = abs(min(x)) + 1;
		ret$y = x + ret$shift;
	}
	return(ret);
}			
			
adjustTrait = function(t, covar, thr = 0.01, marg.thr = 0.01, id.fname='external_name'){
	if(!is.data.frame(covar)) stop("Covariates object must be a data frame");
	if(!is.data.frame(t)) stop("Trait object must be a data frame");
	if(!all(c(id.fname %in% names(covar), id.fname %in% names(t)))) stop(paste("Both input data frames must have an identifier field named:", id.fname));
	c.nms = setdiff(names(covar),id.fname);
	y.nms = unique(t$trait_name);
	if(length(y.nms) != 1) stop("can only adjust for a single trait at a time");
	t = t[t[,id.fname] %in% covar[,id.fname],];
	d = merge(t,covar);
	#ret = t[0,];
	ret = data.frame(matrix(nrow=0,ncol=6));
	names(ret) = c(id.fname,'trait_name','trait_value','shift','adjustment','transform');
	
	#compute marginal pvalues
	shift.val = shiftTrait(d$trait_value);
	d$trait_value = shift.val$y;
	d$shift = shift.val$shift;
	p = c();
	lambda = c();
	for(n in c.nms){
		model = formula(paste(c("trait_value",n),collapse="~"));
		fit = lm(model, data=d,x=T,y=T,model=T);
		bxcx = boxcox(fit,plotit=F);
		p = c(p,anova(fit)[["Pr(>F)"]][1]);
		lambda = c(lambda,bxcx$x[which.max(bxcx$y)]);
	}
	
	#if there are no significant marginal effects, tranform only
	mk = p < marg.thr;
	if(sum(mk) == 0){
		for(l in unique(getLambda(lambda)[getLambda(lambda) != 1])){
			tmp = d[names(t)];
			shift.val = shiftTrait(tmp$trait_value);
			tmp$trait_value = shift.val$y;
			tmp$shift = shift.val$shift;
			tmp$trait_value = bxcxTrans(tmp$trait_value,l);
			tmp$transform = getTransName(l);
			tmp$shift = shift.val$shift;
			tmp$adjustment = 'None';
			
			ret = rbind(ret,tmp);
		}
		
	}else{
		#find the best model
		adjust = paste(sort(c.nms[mk]),collapse='+',sep='');
		model = formula(paste(c("trait_value",adjust),collapse="~"));
		fit = lm(model,data=d,x=T,y=T,model=T);
		fit = stepAIC(fit,trace=0);
		adjust = paste(attr(fit$terms,'term.labels'),collapse='+');
		if(adjust == '') {adjust = 'None';}
		
		#throw out terms with pvalues > threshold
		ind = which(anova(fit)[["Pr(>F)"]] > thr);
		ind1 = which(anova(fit)[["Pr(>F)"]] < thr);
		if(length(ind) > 0 & length(ind1) > 0){
			adjust = paste(sort(attr(fit$terms[-ind],'term.labels')),collapse='+',sep='');
			model = formula(paste(c("trait_value",adjust),collapse="~"));
			fit = update(fit,model);
		}
		bxcx = boxcox(fit,plotit=F);
		
		#store the adjusted traits if adjustment is necessary
		if(adjust != 'None'){
			tmp = d[names(t)];
			shift.val = shiftTrait(tmp$trait_value);
			tmp$trait_value = shift.val$y;
			tmp$shift = shift.val$shift;
			tmp$trait_value = array(fit$resid) +  as.numeric(fit$coefficients['(Intercept)']);
			tmp$shift = shift.val$shift;
			tmp$adjustment = adjust;
			if(!any(names(tmp) == 'transform')){
				tmp$transform = 'normal';
			}
			ret = rbind(ret,tmp);
		}
		
		l= getLambda(bxcx$x[bxcx$y == max(bxcx$y)]);
		if(l != 1){
			d$trait_value = bxcxTrans(d$trait_value,l);
			d$transform = getTransName(l);
			tmp = d[names(t)];
			if(!('shift' %in% names(tmp))) tmp$shift = 0;
			#store the transformed trait prior to adjustment
			ret = rbind(ret,tmp);
		}
		if(l != 1  & length(ind1) > 0){
			adjust = paste(sort(attr(fit$terms,'term.labels')),collapse='+',sep='');
			model = formula(paste(c("trait_value",adjust),collapse="~"));
			fit = lm(model,data=d,x=T,y=T,model=T);	
			fit = stepAIC(fit,trace=0);
			
			#throw out terms with pvalues > threshold
			ind = which(anova(fit)[["Pr(>F)"]] >= thr);
			ind1 = which(anova(fit)[["Pr(>F)"]] < thr);
			if(length(ind) > 0 & length(ind1) > 0){
				adjust = paste(attr(fit$terms[-ind],'term.labels'),collapse='+');
				model = formula(paste(c("trait_value",adjust),collapse="~"));
				fit = update(fit,model);
			}
			pvals = anova(fit)[["Pr(>F)"]];
			ind = which(pvals < thr);
			
			#if significant terms are found, adjust and save the result.
			if(length(ind) > 0){
				adjust = paste(sort(attr(fit$terms,'term.labels')),collapse='+',sep='');
				tmp$trait_value = array(fit$resid) +  as.numeric(fit$coefficients['(Intercept)']);
				tmp$adjustment = adjust;
				ret = rbind(ret,tmp);
			}
		}
	}
	return(ret);
}
