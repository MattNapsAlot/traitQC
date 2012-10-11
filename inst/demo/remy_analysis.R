require(MASS);


####
#config variables
####
root = '//FILER1/GENETICS/External_Data_Repository/NU_Turek/Remy/data';
#root = 'c:/Remy/data';
setwd(root);
source(paste(c(root,'utility.R'),collapse='/'));
fname.trait = paste(c(root,"N2_Remy_phenotypes7_10_08.csv"),collapse='/');
fname.covar = paste(c(root,"RemyCov1.txt"),collapse='/');
pval.thr = 0.01;
id.fname='id';
adj.groups = c('timepoint');

####
#end config
####

#load traits
trait = getTraits(fname.trait,tp.fname='interval',sep=',');
trait = pivotTable(trait,c('external_name','timepoint'));
trait$shift = 0;

#######
#need to clean up timepoint names and compute averages for Sleep studies
#######
# only consider flp, fdp and fld timepoints

trait = trait[!is.na(trait$trait_value),]
trait$timepoint = as.character(trait$timepoint);
ind = c();
for(p in c('flp', 'fld', 'fdp')){
	ind = c(ind, grep(p, trait$timepoint));
}
trait = trait[ind,];

# average metrics
trait$timepoint = sub("[\\.][0-9]$","", trait$timepoint,perl=T);
trait$trait_value = ave(trait$trait_value,trait[,c("external_name","timepoint","trait_name")]);
trait = unique(trait);
rownames(trait) = 1:dim(trait)[1];
trait$adjustment= 'None';
trait$transform='normal';
rm(p, ind);

####
#compute new metrics that are the difference between the light and dark periods
####
tmp = trait[0,];
for(name in unique(trait$external_name)){
	for(tr in unique(trait$trait_name)){
		mk = trait$external_name == name & trait$trait_name == tr & trait$timepoint %in% c('flp','fdp');
		if(sum(mk) != 2 | !all(c('flp','fdp') %in% trait$timepoint[mk])){
			print(paste('Could not compute new metric. Missing data for id:',name,'trait:',tr,collapse=" "));
			next;
		}
		row = dim(tmp)[1] + 1;
		tmp[row,'external_name'] = name;
		tmp[row,'transform'] = 'normal';
		tmp[row,'adjustment'] = 'None';
		tmp[row,'timepoint'] = 'combined'; #maybe not the best timepoint?
		tmp[row,'trait_name'] = paste(c(tr,'flp','minus','fld'),collapse='_');
		tmp[row,'shift'] = 0;
		flp.val = trait[mk,'trait_value'][trait$timepoint[mk] == 'flp'];
		fdp.val = trait[mk,'trait_value'][trait$timepoint[mk] == 'fdp'];
		tmp[row,'trait_value'] = flp.val - fdp.val;
		
	}
}

#append new traits to the list
trait = rbind(trait,tmp);
rm(mk,tmp,flp.val,fdp.val,row,tr,name);

#######
#end Sleep studies cleanup
#######

#load covariates, and relevel
covar = getCovariates(fname.covar,id.fname=id.fname,relevel=T);
new.trait = trait[0,]
for(tr in unique(trait$trait_name)){
	for(tp in unique(trait$timepoint[trait$trait_name==tr])){
		mk = trait$trait_name == tr & trait$timepoint == tp & !is.na(trait$trait_value);
		tmp = trait[mk,];
		if(dim(tmp)[1] != length(unique(tmp$external_name))) {
			print(paste(c("found multiple records for a single individual",tr,tp),collapse=':'));
			next;
		}
		new.trait = rbind(new.trait, adjustTrait(tmp, covar));# use the default pvalue thresholds
		#trait = rbind(trait,new.trait);
	}
}
trait = rbind(trait,new.trait);
rm(tr,tp,mk,tmp,new.trait);

#get tids and merge with trait data
external_name = unique(trait$external_name);
tid = getTidFromExtId(external_name);
id.map = data.frame(tid,external_name);
id.map = id.map[!is.na(id.map$tid),];

trait = merge(id.map,trait);
write.table('all_traits.txt',sep='\t',quote=F,rownames=F);



