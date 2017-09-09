

library(mvtnorm)
library(Matrix)

expit <- function(x){
	exp(x)/(1+exp(x))
	}

gen_null <- function(eg,mea,samp_size=1000){
	#un_norm_vecs <- eg$vectors%*%diag(eg$value)
	un_norm_vecs <- eg$vectors
	#deltas2 <- (t(t(ef_size*eg$value)%*%eg$vectors)^2)
	#deltas2 <- c(t(t(mea)%*%eg$vectors)^2)
	#deltas2 <- (t(t(ef_size)%*%t(un_norm_vecs))^2)
	
	#sqrt_inv_cov <- sqrtm(solve(cov_mat2))
	
	#deltas2 <- (un_norm_vecs%*%sqrt_inv_cov%*%ef_size)^2
	deltas2 <- ((t(un_norm_vecs)%*%mea)/sqrt(eg$values))^2
	
	#tmp <- rep(0,samp_size)
		tmp <- list()
		for (i in seq_along(eg$value)){
			tmp[[i]] <- (eg$values[i])*rchisq(samp_size,df=1,ncp=deltas2[i])
		}
	#print(deltas2)
	vec <- rowSums(matrix(unlist(tmp),nrow=samp_size))
	return(vec)
	}


calc_emp_p <- function(dist,num){
	L <- length(dist)	
	plac <- which(order(c(dist,num))==(L+1))
	pval <- 1-plac/(L+2)
	return(pval)	
	}





cov_mat <- bdiag(list(matrix(0.7,nrow=20,ncol=20) + diag(0.3,20),diag(1,980)))

#table(as.vector(cov_mat[1:500,991:1000]))
# samples along columns per expression norm
prof <- t(rmvnorm(280,rep(8,1000),as.matrix(cov_mat)))
test_mat2 <- cbind(c(rep(1,20),rep(0,980)),c(rep(0,980),rep(1,20)))
colnames(test_mat2) <- c('cor','uncor')

test_list <- list('cor' = 1:20,'uncor' = 981:1000,'test_syntax' = 11:30)


res_lis <- list()

mixed_lis <- list()
dir_lis <- list()

for(i in 1:100){

j <- floor((i-1)/20)
#eff_size <- exp(j/30)-exp(1/30)
eff_size <- 0.3
#eff_size <- 0
lin_pred_uncor <- eff_size*(colSums(prof[990:991,]) -16 )
lin_pred_cor <- eff_size*(colSums(prof[10:11,])  - 16 )

cancer_cor <- rbinom(280,1,prob=expit(lin_pred_cor))
cancer_uncor <- rbinom(280,1,prob=expit(lin_pred_uncor))

##roast_cor <- limma::roast(y= prof, index=1:100,design=cbind(rep(1,length(cancer_cor)),cancer_cor))
roast_cor <- limma::roast(y= prof, index=test_list,design=cbind(rep(1,length(cancer_cor)),cancer_cor),contrast=2,sort='none',nrot=5000) # NOT SURE USING INDEX VECTOR CORRECTLY 

##roast_uncor <- limma::roast(y= prof, index=101:200,design=cbind(rep(1,length(cancer_uncor)),cancer_uncor))
roast_uncor <- limma::roast(y= prof, index=test_list,design=cbind(rep(1,length(cancer_uncor)),cancer_uncor),contrast=2,sort='none',nrot=5000) # NOT SURE USING INDEX VECTOR CORRECTLY 

#data2 <- data[,-c(1,2),with=F]
design_cor <- model.matrix(~ cancer_cor)
fit_cor <- lmFit(prof, design=design_cor)
fit_cor <- eBayes(fit_cor)
#tt <- topTable(fit, coef=2)
#tt

#data2 <- data[,-c(1,2),with=F]
design_uncor <- model.matrix(~ cancer_uncor)
fit_uncor <- lmFit(prof, design=design_uncor)
fit_uncor <- eBayes(fit_uncor)

test_stat <- sum(fit_cor$t[test_list[['cor']],2])
test_stat_uncor <- sum(fit_cor$t[test_list[['uncor']],2])

test_uncor_stat <- sum(fit_uncor$t[test_list[['cor']],2])
test_uncor_stat_uncor <- sum(fit_uncor$t[test_list[['uncor']],2])

test_stat2 <- sum((fit_cor$t[test_list[['cor']],2])^2)
test_stat_uncor2 <- sum((fit_cor$t[test_list[['uncor']],2])^2)

test_uncor_stat2 <- sum((fit_uncor$t[test_list[['cor']],2])^2)
test_uncor_stat_uncor2 <- sum((fit_uncor$t[test_list[['uncor']],2])^2)

look_cor <- cor(t(prof[test_list[['cor']],]))
look_uncor <- cor(t(prof[test_list[['uncor']],]))

emp_dir1 <- 2*pnorm(-abs(test_stat),mean=0,sd=sqrt(sum(look_cor)))
emp_dir2 <- 2*pnorm(-abs(test_uncor_stat_uncor),mean=0,sd=sqrt(sum(look_uncor)))

2*pnorm(-abs(test_uncor_stat),mean=0,sd=sqrt(sum(look_cor)))
2*pnorm(-abs(test_uncor_stat_uncor),mean=0,sd=sqrt(sum(look_uncor)))

#roast_cor[1:2,c(5,7)]
#roast_cor[1:2,5]
#roast_cor[1:2,7]

### for Mixed test
eig_cor <- eigen(look_cor)
eig_uncor <- eigen(look_uncor)

nul_dis_cor <- gen_null(eig_cor,rep(0,dim(look_cor)[1]),samp_size=5000)
nul_dis_uncor <- gen_null(eig_uncor,rep(0,dim(look_uncor)[1]),samp_size=5000)

emp1 <- calc_emp_p(nul_dis_cor,test_stat2)
emp2 <- calc_emp_p(nul_dis_uncor,test_uncor_stat_uncor2)
##print(cbind(roast_cor[1:2,5],c(emp_dir1,emp_dir2)))

# print(cbind(c(roast_cor[1:2,5][1],roast_uncor[1:2,5][2]),c(emp_dir1,emp_dir2)))
# print(cbind(c(roast_cor[1:2,7][1],roast_uncor[1:2,7][2]),c(emp1,emp2)))

dir_lis[[i]] <- cbind(c(roast_cor[1:2,5][1],roast_uncor[1:2,5][2]),c(emp_dir1,emp_dir2))
mixed_lis[[i]] <- cbind(c(roast_cor[1:2,7][1],roast_uncor[1:2,7][2]),c(emp1,emp2))

print(i)
}


dir_cols <- do.call(rbind,dir_lis)
mix_cols <- do.call(rbind,mixed_lis)

cor(mix_cols)
cor(dir_cols)

cor(-log(mix_cols))
cor(-log(dir_cols))

pdf('~/Google_Drive/Oslo_sync/gsea_paper/alt_mix.pdf')
plot(mix_cols[,1],mix_cols[,2],cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="p-val Roast",ylab="p-val Chi-sq mixture",main='Mixed Test')
abline(0,1)
dev.off()

pdf('~/Google_Drive/Oslo_sync/gsea_paper/alt_dir2.pdf')
plot(dir_cols[,1],dir_cols[,2],cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="p-val Roast",ylab="p-val Chi-sq mixture")
# ,main='Directed Test'
abline(0,1)
dev.off()

pdf('~/Google_Drive/Oslo_sync/gsea_paper/alt_mix_log.pdf')
plot(-log(mix_cols[,1]),-log(mix_cols[,2]),cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="-log p Roast",ylab="-log p Chi-sq mixture",main='Mixed Test')
abline(0,1)
dev.off()

pdf('~/Google_Drive/Oslo_sync/gsea_paper/alt_dir_log.pdf')
plot(-log(dir_cols[,1]),-log(dir_cols[,2]),cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="-log p Roast",ylab="-log p Chi-sq mixture",main='Directed Test')
abline(0,1)
dev.off()


save(dir_cols,mix_cols,file='~/Google_Drive/Oslo_sync/gsea_paper/alt_roast_similarity.RObject',compress=F)


load('~/Google_Drive/Oslo_sync/gsea_paper/alt_roast_similarity.RObject')
colMeans(dir_cols)
colMeans(mix_cols)

load('~/Google_Drive/Oslo_sync/gsea_paper/null_roast_similarity.RObject')
colMeans(dir_cols)
colMeans(mix_cols)



pdf('~/Google_Drive/Oslo_sync/gsea_paper/null_mix.pdf')
plot(mix_cols[,1],mix_cols[,2],cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="p-val Roast",ylab="p-val Chi-sq mixture",main='Null Mixed Test')
abline(0,1)
dev.off()

pdf('~/Google_Drive/Oslo_sync/gsea_paper/null_dir2.pdf')
plot(dir_cols[,1],dir_cols[,2],cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="p-val Roast",ylab="p-val Chi-sq mixture")
#,main='Null Directed Test'
abline(0,1)
dev.off()

pdf('~/Google_Drive/Oslo_sync/gsea_paper/null_mix_log.pdf')
plot(-log(mix_cols[,1]),-log(mix_cols[,2]),cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="-log p Roast",ylab="-log p Chi-sq mixture",main='Null Mixed Test')
abline(0,1)
dev.off()

pdf('~/Google_Drive/Oslo_sync/gsea_paper/null_dir_log.pdf')
plot(-log(dir_cols[,1]),-log(dir_cols[,2]),cex.lab=1.5,cex.axis=1.5,cex.main=1.45,xlab="-log p Roast",ylab="-log p Chi-sq mixture",main='Null Directed Test')
abline(0,1)
dev.off()

save(dir_cols,mix_cols,file='~/Google_Drive/Oslo_sync/gsea_paper/null_roast_similarity.RObject',compress=F)










cov_mat2 <- cov_mat[1:40,1:40]
ran_me <- abs(rnorm(40))

prof <- rmvnorm(1000,mean=ran_me,sigma = as.matrix(cov_mat2))


sv <- svd(prof)
delta <- rep(99,40)

for( i in seq(40))
	delta[i] <- crossprod(sv$v[,i],ran_me) / (sv$d[i])


samp <- matrix(-99,nrow=1000,ncol=40)
for(i in seq(40))
	samp[,i] <- sv$d[i]*rchisq(1000,df=1,ncp=delta[i]^2)


cor_samps <- rowSums(prof^2)
sim_samps <- rowSums(samp)


qqplot(cor_samps,sim_samps)


eig_exam <- eigen(cov_mat2)
#est_cov <- crossprod(prof,prof)/(dim(prof)[1])
est_cov <- cov(prof)
eig_exam <- eigen(est_cov)

sim_samp2 <- gen_null(eig_exam,ran_me)

qqplot(cor_samps,sim_samp2)
abline(0,1)



