# This script uses two-way t-test to infer regulatory relationships between candidate regulators and genes.
# Input 1. a matrix of expression values for genes under a  set of conditions. Each row is the expression values of a gene,
#	  while each column denotes the expression values of all genes under a particular condtions. 
#	  Hence, row names are genes and column names are samples. Expression values are 
#	  seperaed by tabs(\t). 
#      2. a list of candidate transcription factors whose expression values are included in the frist parameter
#      3. output file name
# Output 1. a ranked list of predicted reg-gene interactions
#
# Copyright (C) 2012, Jianlong Qi and Tom Michoel
#
#   Jianlong Qi <jianlong.qi@frias.uni-freiburg.de>
#   Tom Michoel <tom.michoel@roslin.ed.ac.uk>
#
#   http://www.roslin.ed.ac.uk/tom-michoel
#


rm(list=ls());
library(limma);
#input parameters
expression.file = "all.expression.normalized.tab"; 
regulator.file = "regulator.tab"; 
output.file ="regulatory_relationship.tab"

#read expression values and regulators
expression = read.table(file=expression.file,sep="\t",quote="",comment.char="");
regulators = read.table(file=regulator.file,sep="\t",
			quote="",comment.char="",stringsAsFactors=F)[,1];


#Determine the contrasts where genes are most differentially expressed 
#using moderated t-statistics
order.expression = apply(expression,1,order);
sort.expression = t(apply(expression,1,sort,na.last=T));
ordinary.t.test = numeric();
#Loop for calculating the moderated t-statistics between 
#the top i and remain conditions for each gene
for (i in 1:(dim(sort.expression)[2]-1)){
	t.con.clusters = rep(1,dim(sort.expression)[2]);
        t.con.clusters[1:i] = 2;
        design = model.matrix(~0+factor(t.con.clusters));
        colnames(design) = paste("c",seq(1,2),sep="");
        contrast.matrix = makeContrasts(c1-c2,levels=design);
        fit = lmFit(sort.expression,design);
        fit2 = contrasts.fit(fit,contrast.matrix);
        fit3 = eBayes(fit2);
        moderated.t.test = abs(fit3$t[,1]); 
        moderated.t.test[which(is.na(fit3$t[,1])==T)] = 0;
	ordinary.t.test = cbind(ordinary.t.test,moderated.t.test);
}

#Select differentially expressed transcription factors in the critical contrast
#of each gene. The ranked list of TFs are saved as the output.
gene.critical.contrast = apply(t(ordinary.t.test),2,which.max);
regulator.expression = expression[regulators,];
reglatory.relation = data.frame();#data frame for regulatory relationships
for ( i in 1:length(gene.critical.contrast)){
	#the gene for the current iteration
	gene.name = names(gene.critical.contrast)[i];
	#build contrast matrix for the critical contrast
	critical.contrast.pivot = gene.critical.contrast[gene.name];
	critical.contrast = rep(1,dim(regulator.expression)[2]);
        critical.contrast[order.expression[1:critical.contrast.pivot,gene.name]] = 2;
        design = model.matrix(~0+factor(critical.contrast));
        colnames(design) = paste("c",seq(1,2),sep="");
        contrast.matrix = makeContrasts(c1-c2,levels=design);
	#apply LIMMA to find differentailly expressed transcription factors in the contrast
        fit = lmFit(regulator.expression,design);
        fit2 = contrasts.fit(fit,contrast.matrix);
        fit3 = eBayes(fit2);
        moderated.t.test = abs(fit3$t[,1]); #use t-stat
        moderated.t.test[which(is.na(fit3$t[,1])==T)] = 0;
	#normalized t-statistics
	reg.score = (moderated.t.test-mean(moderated.t.test,na.rm=T))/sd(moderated.t.test,na.rm=T);
	reglatory.relation = rbind(reglatory.relation,
	                           data.frame(regulator=names(reg.score),
				   gene=gene.name,score=reg.score[names(reg.score)]));	
}

#save result
write.table(reglatory.relation,file=output.file,sep="\t",quote=F);	
