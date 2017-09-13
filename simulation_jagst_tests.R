


library(mvtnorm)
library(Matrix)
library(selectiveInference)


expit <- function(x){
	exp(x)/(1+exp(x))
	}

brute_force_null <- function(stats,cov,siz,its){
		len <- length(stats)
		null <- rep(-99,its)
	for(kk in seq(its)){
		inds <- sample(1:len,siz,replace=F)
		null[kk] <- t(stats[inds])%*%solve(cov[inds,inds])%*%(stats[inds])
		}
		return(null)
	}

calc_emp_p <- function(dist,num){
	L <- length(dist)	
	plac <- which(order(c(dist,num))==(L+1))
	pval <- 1-plac/(L+2)
	return(pval)	
	}

### MIXTURE OF CHI-SQUARES WITH DIFFERENT NCPS? -- use HARVILLE PAPER

easier_null <- function(len_inds,dim_x,its,num_nest,cancer_cor){
	stat_nulls <- list()
	for(ee in 1:its){
	null_sam_inds <- sample(1:dim_x,len_inds,replace=F)
	
	larfit = lar(mat_trans[,null_sam_inds],cancer_cor) # is there a way to stop the algorithm early in general since aic will likely be achieve early in the calculation
	  out.aic <- larInf(larfit,type="aic")
	  ncps_cor <- round(out.aic$sign * out.aic$vmat %*% out.aic$y/(out.aic$sigma * sqrt(rowSums(out.aic$vmat^2))),3) # test stats
		# test_stat_cor,
		
	  stat_nulls[[ee]] <- rchisq(num_nest,df=length(inds1),ncp=sum(ncps_cor^2)) ## in these simulations would generally be assuming same size of correlated and uncorrelated test set so should not necessarily need to change df parameter but doing so here to prevent future issues
		}
		return(unlist(stat_nulls))	
	}





mat_trans <- t(as.matrix(prof))
len_inds <- 20
dim_x <- dim(mat_trans)[2]
its <- 25
num_nest <- 200

cov_mat <- bdiag(list(matrix(0.7,nrow=20,ncol=20) + diag(0.3,20),diag(1,20)))

#table(as.vector(cov_mat[1:500,991:1000]))
# samples along columns per expression norm
prof <- t(rmvnorm(280,rep(8,40),as.matrix(cov_mat)))
test_mat2 <- cbind(c(rep(1,20),rep(0,20)),c(rep(0,20),rep(1,20)))
colnames(test_mat2) <- c('cor','uncor')

test_list <- list('cor' = 1:20,'uncor' = 21:40)
res_lis <- list()

inds1 <- test_list[[1]] 
inds2 <- test_list[[2]]

#for(i in 1:80){
for(i in 1:200){
	#j <- floor((i-1)/4)

	# divide by 3 for self-contained test
	j <- floor((i-1)/10)/2

	eff_size <- exp(j/30)-exp(1/30)
	#eff_size <- 0
	lin_pred_uncor <- eff_size*(colSums(prof[30:31,]) - 16 )
	lin_pred_cor <- eff_size*(colSums(prof[10:11,])  - 16 )
	
	cancer_cor <- rbinom(280,1,prob=expit(lin_pred_cor))
	cancer_uncor <- rbinom(280,1,prob=expit(lin_pred_uncor))
	
	des <- cbind(1,cancer_cor)
	lmOb <- limma::lmFit(prof,design=des)
	lmOb2 <- limma::eBayes(lmOb)
	statsCor <- lmOb2$t[,2]
	# unmoderated for self-contained test
	  statsCor <- (lmOb2$coef / lmOb2$stdev.unscaled / lmOb2$sigma)[,2]


	des <- cbind(1,cancer_uncor)
	lmOb <- limma::lmFit(prof,design=des)
	lmOb2 <- limma::eBayes(lmOb)
	statsUncor <- lmOb2$t[,2]
	# unmoderated for self-contained test
	  statsUncor <- (lmOb2$coef / lmOb2$stdev.unscaled / lmOb2$sigma)[,2]

	
	test_stat_cor <- t(statsCor[inds1])%*%solve(cov_mat[inds1,inds1])%*%(statsCor[inds1])
	test_stat_uncor <- t(statsUncor[inds2])%*%solve(cov_mat[inds2,inds2])%*%(statsUncor[inds2])
	
	# null_distCor <- brute_force_null(statsCor,cov_mat,20,1000)
	# null_distUncor <- brute_force_null(statsUncor,cov_mat,20,1000)

	# brute_cor <- calc_emp_p(null_distCor,test_stat_cor[1])
	# brute_uncor <- calc_emp_p(null_distUncor,test_stat_uncor[1])
	
	#### NEW CALC #####
	te_null_cor <- easier_null(len_inds,dim_x,its,num_nest,cancer_cor)
	te_null_uncor <- easier_null(len_inds,dim_x,its,num_nest,cancer_uncor)

	p_self_cor <- pchisq(test_stat_cor[1],df=(len_inds),lower.tail=F)
	p_self_uncor <- pchisq(test_stat_uncor[1],df=(len_inds),lower.tail=F)

	 # larfit = lar(x,y) # is there a way to stop the algorithm early in general since aic will likely be achieve early in the calculation
		# out.aic = larInf(larfit,type="aic")
		# x <- out.aic
		# ncps_cor <- round(x$sign * x$vmat %*% x$y/(x$sigma * sqrt(rowSums(x$vmat^2))),3) # test stats
			
		# p <- pchisq(test_stat_cor,df=length(inds1),ncp=sum(ncps_cor^2),lower.tail=F) ## in these simulations would generally be assuming same size of correlated and uncorrelated test set so should not necessarily need to change df parameter but doing so here to prevent future issues
		
		# p <- pchisq(test_stat_uncor,df=length(inds2),ncp=sum(ncps_uncor^2),lower.tail=F) 		
	easy_cor <- calc_emp_p(te_null_cor,test_stat_cor[1])
	easy_uncor <- calc_emp_p(te_null_uncor,test_stat_uncor[1])

	#### END OF NEW CALC #####

	#safe::safe(as.matrix(prof),cancer_cor,asdfasdf)
	
	# (stats,cov,siz,its){
	# brute_cor <- brute_force_null(stats,cov_mat,cancer_cor,1000)
	# brute_uncor <- brute_force_null(stats,cov_mat,cancer_uncor,1000)
	
	# safe_cor@global.pval
	# safe_uncor@global.pval
	
	# eleg_cor <- safe::safe(as.matrix(prof),cancer_cor,asdfasdf)
	# eleg_uncor <- safe::safe(as.matrix(prof),cancer_uncor,sdfgsdfg0)
	# tmp_res <- list('bruteCor'= roast_cor,'bruteUncor'= roast_uncor,'elegCor'=camera_cor,'elegUncor'=camera_uncor)
	
	##tmp_res <- list('bruteCor'= brute_cor,'bruteUncor'= brute_uncor,'easyCor'= easy_cor, 'easy_uncor' = easy_uncor)
	tmp_res <- list('easyCor'= easy_cor, 'easyUncor' = easy_uncor,'self_p_cor' = p_self_cor,'self_p_uncor' = p_self_uncor)
	
	res_lis[[i]] <- tmp_res
	print(tmp_res)
	print(i)
}


##no_force_null <- function(stats,cov,siz,its){
	






plot(lowess(seq(20),colMeans(matrix(-log(sapply(res_lis,function(x){x$self_p_cor})),nrow=10)),f=2/3),col='black',ty='l',ylim=c(0,8),lwd=4,cex.axis=1.4,cex.lab=1.8,xlab="effect size of two transcripts",ylab='-log p-value')
lines(lowess(seq(20),colMeans(matrix(-log(sapply(res_lis,function(x){x$self_p_uncor})),nrow=10)),f=2/3),col='red',lwd=4)


plot(lowess(seq(20),colMeans(matrix(-log(1-sapply(res_lis,function(x){x$easyCor})),nrow=4)),f=2/3),col='black',ty='l',ylim=c(0,8),lwd=4,cex.axis=1.4,cex.lab=1.8,xlab="effect size of two transcripts",ylab='-log p-value')
lines(lowess(seq(20),colMeans(matrix(-log(1-sapply(res_lis,function(x){x$easy_uncor})),nrow=4)),f=2/3),col='red',lwd=4)

