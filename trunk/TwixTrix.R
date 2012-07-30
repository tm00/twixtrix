twixtrix = function (exp, tfs = NULL){
# TWIXTRIX - Double 2-way t-test network inference
# TWIXTRIX implements the double 2-way t-test network inference method
# using unmoderated t-statistics
#
# INPUT:
#   - exp: data matrix, rows are genes, columns are samples
#   - tfs: list of regulators in the form of row numbers
#
# OUTPUT:
#   - score : matrix of size #tfs x #genes containing t-statistic for each
#             TF in the critical contrast of each gene
#   - Z     : Z-transcformed score
# Copyright (C) 2012, Jianlong Qi and Tom Michoel
#   Jianlong Qi <jianlong.qi@frias.uni-freiburg.de>
#   Tom Michoel <tom.michoel@roslin.ed.ac.uk>
#   http://www.roslin.ed.ac.uk/tom-michoel
	
	K = dim(exp)[1] #number of genes
	N = dim(exp)[2] #number of samples
	
	x = t(apply(exp, 1, sort))	#sort data rows
	x2 = x * x
	c_x = t(apply(x, 1, cumsum))	#cumulative row sums
	c_x2 = t(apply(x2, 1,cumsum))	#cumulative sum of squares
	
	#mean & std for group 1 and 2;
	#note that s = (n-1)*sd^2 for group 1 and 2
	n1 = matrix(rep(seq(1, N-1), times = K), nrow = K, byrow = T)
	n2 = N - n1
	u1 = c_x[, -N] / n1
	u2 = (c_x[, N] - c_x[, -N]) / n2
	s1 = c_x2[, -N] - n1 * u1 * u1
	s2 = c_x2[, N] - c_x2[, -N] - n2 * u2 * u2

	#all t-statistics
	ot = abs(u1 - u2) / sqrt((s1 + s2) * N / ((N-2) * n1 * n2))
	
	
	#critical contrast
	o = t(apply(exp, 1, order))
	max_ot = apply(ot, 1, which.max)
	p = matrix(F, nrow = K, ncol = N)
	p = mapply(function(a,b,c){a[b[1:c]]=T; return(a)}, data.frame(t(p)), data.frame(t(o)), max_ot)
	
	if (is.null(tfs)) tfs = seq(1, K) #all against all  
	exp_tfs = exp[tfs, ]
	exp_tfs2 = exp_tfs * exp_tfs
	
	n1_tf = colSums(p)
	u1_tf = t(exp_tfs %*% p) / n1_tf 
	s1_tf = t(exp_tfs2 %*% p) - n1_tf * u1_tf * u1_tf
	
	n2_tf = colSums(1-p)
	u2_tf = t(exp_tfs %*% (1-p)) / n2_tf
	s2_tf = t(exp_tfs2 %*% (1-p)) - n2_tf * u2_tf * u2_tf
	
	#t-statistics of TFs in critical contrast
	ot_tf = abs(u1_tf - u2_tf) / sqrt((s1_tf + s2_tf) * N / ((N-2) * n1_tf * n2_tf))
	
	#remove self-interactions
	for (tf  in 1:length(tfs)) {ot_tf[tfs[tf], tf] = 0}
	
	#normalize over genes
	Z_score = scale(t(ot_tf))

	return(list(score=t(ot_tf), Z=Z_score))
}
