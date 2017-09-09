
library(mvtnorm)
library(Matrix)

expit <- function(x){
	exp(x)/(1+exp(x))
	}

cov_mat <- bdiag(list(matrix(0.7,nrow=20,ncol=20) + diag(0.3,20),diag(1,980)))

#table(as.vector(cov_mat[1:500,991:1000]))
# samples along columns per expression norm
prof <- t(rmvnorm(280,rep(8,1000),as.matrix(cov_mat)))
test_mat2 <- cbind(c(rep(1,20),rep(0,980)),c(rep(0,980),rep(1,20)))
colnames(test_mat2) <- c('cor','uncor')

test_list <- list('cor' = 1:20,'uncor' = 981:1000,'test_syntax' = 11:30)


prof_mat <- as.matrix(prof)


gen_names <- paste0(rep(LETTERS,each=50),1:50)[1:1000]
rownames(prof_mat) <- gen_names
colnames(prof_mat) <- 1:280

res_lis <- list()
for(i in 1:400){

j <- floor((i-1)/20)
eff_size <- exp(j/30)-exp(1/30)
#eff_size <- 0
lin_pred_uncor <- eff_size*(colSums(prof[990:991,]) -16 )
lin_pred_cor <- eff_size*(colSums(prof[10:11,])  - 16 )

cancer_cor <- rbinom(280,1,prob=expit(lin_pred_cor))
cancer_uncor <- rbinom(280,1,prob=expit(lin_pred_uncor))

##roast_cor <- limma::roast(y= prof, index=1:100,design=cbind(rep(1,length(cancer_cor)),cancer_cor))
# # # # # # # roast_cor <- limma::roast(y= prof, index=test_list,design=cbind(rep(1,length(cancer_cor)),cancer_cor),contrast=2,sort='none') # NOT SURE USING INDEX VECTOR CORRECTLY 

##roast_uncor <- limma::roast(y= prof, index=101:200,design=cbind(rep(1,length(cancer_uncor)),cancer_uncor))
# # # # # # roast_uncor <- limma::roast(y= prof, index=test_list,design=cbind(rep(1,length(cancer_uncor)),cancer_uncor),contrast=2,sort='none') # NOT SURE USING INDEX VECTOR CORRECTLY 

##camera_cor <- limma::camera(prof, seq(100),design=cbind(rep(1,length(cancer_cor)),cancer_cor))
# # # # # # # # # camera_cor <- limma::camera(prof, test_list,design=cbind(Intercept = rep(1,length(cancer_cor)),Group = cancer_cor),contrast=2,inter.gene.cor=NA,sort=FALSE)
##camera_cor <- limma::camera(prof, index = 1:100,design=cbind(Intercept = rep(1,length(cancer_cor)),Group = cancer_cor),contrast=2)

##camera_uncor <- limma::camera(prof, seq(101,200),design=cbind(rep(1,length(cancer_uncor)),cancer_uncor))

# # # # # # # # # camera_uncor <- limma::camera(prof, test_list ,design=cbind(Intercept = rep(1,length(cancer_uncor)),Group = cancer_uncor),contrast=2,inter.gene.cor=NA,sort=FALSE)

camera_cor_ranks <- limma::camera(prof, test_list ,design=cbind(rep(1,length(cancer_uncor)),cancer_cor),contrast=2,use.ranks=T,inter.gene.cor=NA,sort=FALSE)

camera_uncor_ranks <- limma::camera(prof, test_list ,design=cbind(rep(1,length(cancer_uncor)),cancer_uncor),contrast=2,use.ranks=T,inter.gene.cor=NA,sort=FALSE)

#limma::cameraPR(data2, inds_gdc(inflamm),design=smoke)

# # # # # # # # safe_cor <- safe::safe(as.matrix(prof),cancer_cor,C.mat=test_mat2,Pi.mat=1000)
# # # # # # # # safe_uncor <- safe::safe(as.matrix(prof),cancer_uncor,C.mat=test_mat2,Pi.mat=1000)

# safe_cor@global.pval
# safe_uncor@global.pval



### NEW GAGE WORK  ####
	gage_cor <- gage::gage(prof_mat, gsets = list(cor = gen_names[1:20]),ref = which(cancer_cor==0),same.dir=F,samp = which(cancer_cor==1),use.fold=F,compare = 'unpaired')
#	gage_cor <- gage::gage(gage_mat[,which(typ_blood==type)], gsets = list(first = as.character(data[[2]][which(gsea_list[[gen_ind]])])),ref = which(smoke[which(typ_blood==type)]==0),same.dir=F,samp = which(smoke[which(typ_blood==type)]==1),use.fold=F,compare = 'unpaired')

	gage_uncor <- gage::gage(prof_mat, gsets = list(uncor = gen_names[981:1000]),ref = which(cancer_uncor==0),same.dir=F,samp = which(cancer_uncor==1),use.fold=F,compare = 'unpaired')
#	gage_uncor <- gage::gage(gage_mat[,which(typ_blood==type)], gsets = list(cor = gen_names[1:20],uncor = gen_nams[981:1000]),ref = which(smoke[which(typ_blood==type)]==0),same.dir=F,samp = which(smoke[which(typ_blood==type)]==1),use.fold=F,rank.test=T,compare = 'unpaired')

	gage_rank_cor <- gage::gage(prof_mat, gsets = list(cor = gen_names[1:20]),ref = which(cancer_cor==0),rank.test = T,same.dir=F,samp = which(cancer_cor==1),use.fold=F,compare = 'unpaired')

	gage_rank_uncor <- gage::gage(prof_mat, gsets = list(uncor = gen_names[981:1000]),rank.test = T,ref = which(cancer_uncor==0),same.dir=F,samp = which(cancer_uncor==1),use.fold=F,compare = 'unpaired')
### END OF NEW GAGE WORK  ####


#########tmp_res <- list('roastCor'= roast_cor,'roastUncor'= roast_uncor,'cameraCor'=camera_cor,'cameraUncor'=camera_uncor,'RcameraCor'=camera_cor_ranks,'RcameraUncor'=camera_uncor_ranks,'safeCor'=safe_cor,'safeUncor'=safe_uncor,'gageCor'=gage_cor,'gageUncor'=gage_uncor,'RgageCor'=gage_rank_cor,'RgageUncor'=gage_rank_uncor)
tmp_res <- list('RcameraCor'=camera_cor_ranks,'RcameraUncor'=camera_uncor_ranks,'gageCor'=gage_cor,'gageUncor'=gage_uncor,'RgageCor'=gage_rank_cor,'RgageUncor'=gage_rank_uncor)

res_lis[[i]] <- tmp_res

print(i)
}

#save(res_lis,file='~/Google_Drive/Oslo_sync/power_gsea_methods.RData',compress=FALSE)
save(res_lis, file='~/Google_Drive/Oslo_sync/gsea_paper/power_gsea_gage_camRanks2.RData',compress=FALSE)



fit <- limma::lmFit(prof, design=cbind(rep(1,length(cancer_cor)),cancer_cor),contrast=2)
fit2 <- limma::eBayes(fit)


fit2$p.value

lapply(res_lis,function(x){x$cameraCor})
lapply(res_lis,function(x){x$cameraUncor})
lapply(res_lis,function(x){x$RcameraCor})
lapply(res_lis,function(x){x$RcameraUncor})

lapply(res_lis,function(x){x$roastCor})
lapply(res_lis,function(x){x$roastUncor})

lapply(res_lis,function(x){x$safeCor@global.pval})
lapply(res_lis,function(x){x$safeUncor@global.pval})




limma:lmFit()

str(res_lis[[1]][[6]])

res_lis[[1]]$cameraCor
res_lis[[1]]$cameraUncor


rocam_syntax <- function(res_lis,str,row,col,num_reps=5){
	colMeans(matrix(-log(unlist(lapply(res_lis,function(x){x[[str]][[col]][row]}))),nrow=num_reps))
}



eff_size_graph <- exp((1:20)/30)-exp(1/30)

plot(lowess(eff_size_graph,rocam_syntax(res_lis,'cameraCor',1,4,20),f=2/5,delta=0),ty='l',ylim=c(-1,30),cex.lab=1.5,cex.axis=1.3,xlab="effect size of two transcripts",ylab='-log p-value',lwd=5)
lines(lowess(eff_size_graph,rocam_syntax(res_lis,'cameraUncor',2,4,20),f=2/5,delta=0),col='green',lwd=5)
grid(lwd=2)

grid(7,7,lwd=2)


plot(supsmu(eff_size_graph,rocam_syntax(res_lis,'roastCor',1,5,20),span=0.15,wt=c(rep(10,4),rep(1,16))),ty='l',lwd=5,ylim=c(0,7),cex.axis=1.4,cex.lab=1.8,xlab="effect size of two transcripts",ylab='-log p-value')
lines(supsmu(eff_size_graph,rocam_syntax(res_lis,'roastUncor',2,5,20),span=0.4,wt=c(0.5,0.5,0.5,rep(2,17))),col='green',lwd=5)

plot(supsmu(eff_size_graph,safe_syntax(res_lis,'safeCor',1,20),span=0.3,wt=c(rep(10,4),rep(1,16))),ty='l',lwd=5,ylim=c(0,10),cex.axis=1.4,cex.lab=1.8,xlab="effect size of two transcripts",ylab='-log p-value')
lines(supsmu(eff_size_graph,safe_syntax(res_lis,'safeUncor',2,20),span=0.3,wt=c(rep(10,4),rep(1,16))),col='purple',lwd=5)


plot(supsmu(eff_size_graph,rocam_syntax(res_lis,'roastCor',1,5,20),span=0.15,wt=c(rep(10,4),rep(1,16))),ty='l',lwd=5,ylim=c(0,7),cex.axis=1.4,cex.lab=1.8,xlab="effect size of two transcripts",ylab='-log p-value')
lines(supsmu(eff_size_graph,rocam_syntax(res_lis,'roastUncor',2,5,20),span=0.4,wt=c(0.5,0.5,0.5,rep(2,17))),col='green',lwd=5)








plot(lowess(eff_size_graph,safe_syntax(res_lis,'safeCor',1,20),f=1/2),ty='l',,lwd=5,ylim=c(0,10))
lines(lowess(eff_size_graph,safe_syntax(res_lis,'safeUncor',2,20),f=1/2),col='purple',lwd=5)








library(data.table)
library(ggplot2)

dat <- data.table('-log p-value'=rocam_syntax(res_lis,'roastCor',1,5,20),"effect size of two transcripts"=eff_size_graph)


dat <- data.table('y'=rocam_syntax(res_lis,'roastCor',1,5,20),"x"=eff_size_graph)


dat2 <- data.table('y'=c(rocam_syntax(res_lis,'roastCor',1,5,20),rocam_syntax(res_lis,'roastUncor',2,5,20)),"x"=rep(eff_size_graph,2),'gr'=rep(0:1,each=20))



ggplot(data=dat) + geom_smooth(aes('-log p-value' ~ "effect size of two transcripts")) + theme_bw()


ggplot(dat2, aes(x, qchisq( exp(-y),1,lower.tail=F),group=gr)) +
#  geom_point() +
  geom_smooth(span = 0.7,se=F) + theme_bw()



cars.lo <- loess(dist ~ speed, cars)
predict(cars.lo, data.frame(speed = seq(5, 30, 1)), se = TRUE)


plot(lowess(eff_size_graph,rocam_syntax(res_lis,'roastCor',1,5,20),f=1/2,delta=0.000001),ty='l',lwd=5)
lines(lowess(eff_size_graph,rocam_syntax(res_lis,'roastUncor',2,5,20),f=1/2,delta=0.00001),col='green',lwd=5)





plot(loess(rocam_syntax(res_lis,'roastCor',1,5,20) ~ eff_size_graph,span=0.9),ty='l',lwd=5)
lines(lowess(eff_size_graph,rocam_syntax(res_lis,'roastUncor',2,5,20),f=1/2,delta=0.00001),col='green',lwd=5)



plot(supsmu(rocam_syntax(res_lis,'roastCor',1,5,20), eff_size_graph),ty='l',lwd=5)

plot(lowess(-log(unlist(lapply(res_lis,function(x){x[['roastCor']][[5]][1]}))),rep(eff_size_graph,each=20),f=2/3),ty='l',lwd=5)





plot(eff_size_graph,rocam_syntax(res_lis,'roastCor',1,5,20),ty='l')
lines(eff_size_graph,rocam_syntax(res_lis,'roastUncor',2,5,20),col='green')



plot(eff_size_graph,rocam_syntax(res_lis,'roastCor',1,5,5),ty='l')
lines(eff_size_graph,rocam_syntax(res_lis,'roastUncor',2,5,5),col='green')



safe_syntax <- function(res_lis,str,num,num_reps=5){
	colMeans(matrix(-log(unlist(lapply(res_lis,function(x){x[[str]]@global.pval[num]}))),nrow=num_reps))
}



matrix(-log(unlist(lapply(res_lis,function(x){x[['safeCor']]@global.pval[1]}))),nrow=5)





lines(lowess(eff_size_graph,safe_syntax(res_lis,'safeCor',2,20),f=1/3),col='orange',lwd=3)
lines(lowess(eff_size_graph,safe_syntax(res_lis,'safeUncor',1,20),f=1/3),col='blue',lwd=3)

##plot(seq(20),safe_syntax(res_lis,'safeCor',1))








plot(lowess(seq(20),-log(unlist(lapply(res_lis,function(x){x$safeCor@global.pval[1]})))))
plot(lowess(seq(20),-log(unlist(lapply(res_lis,function(x){x$safeCor@global.pval[2]})))))


plot(lowess(seq(20),-log(unlist(lapply(res_lis,function(x){x$safeUncor@global.pval[1]})))))
plot(lowess(seq(20),-log(unlist(lapply(res_lis,function(x){x$safeUncor@global.pval[2]})))))

plot(seq(20),-log(unlist(lapply(res_lis,function(x){x$safeUncor@global.pval[1]}))))
plot(seq(20),-log(unlist(lapply(res_lis,function(x){x$safeUncor@global.pval[2]}))))


summary(-log(unlist(lapply(res_lis,function(x){x$safeUncor@global.pval[1]}))))
summary(-log(unlist(lapply(res_lis,function(x){x$safeUncor@global.pval[2]}))))




res_lis[[1]]$roastCor[,7]  ## using mixed p-value here, synonymous wtith undirected test
res_lis[[1]]$roastUncor[,7]  ## using mixed p-value here, synonymous wtith undirected test

res_lis[[1]]$safeCor@global.pval
res_lis[[1]]$safeUncor@global.pval

str(res_lis[[1]]$roastCor)


(res_lis[[1]][[6]])@global.pval






cov_mat2 <- matrix(0.6,nrow=50,ncol=50) + diag(0.4,50)
y_add <- rmvnorm(100,mean=rep(0,50),sigma=cov_mat2)


for(i in 1:20){
y <- matrix(rnorm(2000*100),2000,100)

y[1:50,] <- t(y_add)
design <- cbind(Intercept=1,Group=c(rep(0,50),rep(1,50)))

# First set of 20 genes are genuinely differentially expressed
index1 <- 1:50
#y[index1,4:6] <- y[index1,4:6]+1

# Second set of 20 genes are not DE
index2 <- 21:40
 
print(limma::camera(y, index1, design,inter.gene.cor=NA))
print(limma::camera(y, index2, design,inter.gene.cor=NA))

}

