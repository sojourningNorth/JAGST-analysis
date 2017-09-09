
###  NO MORE TCGA IN THIS ONE 

library(mvtnorm)
library(Matrix)
library(selectiveInference)

library(data.table)
library(limma)
library(GEOquery)
library(KEGGREST)
library(GO.db)
library(org.Hs.eg.db)
library(corrplot)

library(safe)
library(org.Hs.eg.db)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)

# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
# biocLite("AnnotationHub")


get_entrez_from_kegg <- function(p53_path){
	len_gen <- length(p53_path[[1]]$GENE)
	tmp_ret <- p53_path[[1]]$GENE[(2*seq(l/2)-1)]
	return(tmp_ret)
	}
	
kegg_entrez <- function(p53_path){
	l <- length(p53_path[[1]]$GENE)
	entrez <- p53_path[[1]]$GENE[2*seq(l/2)-1]
	return(entrez)
}

perc_in_gds <- function(inflamm){
	return(mean(unique(inflamm)%in%unlist(sym_eg[unique(as.character(data[[2]]))])))
}

inds_gdc <- function(inflamm,log=F){	
	return(dat2_vals_fin%in%inflamm)
}


expit <- function(x){
	exp(x)/(1+exp(x))
	}
	
calc_emp_p <- function(dist,num){
	L <- length(dist)	
	plac <- which(order(c(dist,num))==(L+1))
	pval <- 1-plac/(L+2)
	return(pval)	
	}

### MIXTURE OF CHI-SQUARES WITH DIFFERENT NCPS? -- use HARVILLE PAPER


easier_null <- function(mat_trans,len_inds,its=200,num_nest=1000,cancer_cor){ 
	stat_nulls <- list()
	
	dim_x <- dim(mat_trans)[2]
	
	for(ee in 1:its){
	null_sam_inds <- sample(1:dim_x,len_inds,replace=F)
	
	larfit = lar(mat_trans[,null_sam_inds],cancer_cor) # is there a way to stop the algorithm early in general since aic will likely be achieve early in the calculation
	  out.aic <- larInf(larfit,type="aic")
	  ncps_cor <- round(out.aic$sign * out.aic$vmat %*% out.aic$y/(out.aic$sigma * sqrt(rowSums(out.aic$vmat^2))),3) # test stats
		# test_stat_cor,
	  stat_nulls[[ee]] <- rchisq(num_nest,df=len_inds,ncp=sum(ncps_cor^2)) ## in these simulations would generally be assuming same size of correlated and uncorrelated test set so should not necessarily need to change df parameter but doing so here to prevent future issues
		}
		return(unlist(stat_nulls))	
	}

polymorph_test <- function(dat_test,out_test,inds1){ # dat_test is a matrix, out_test a vector, inds1 a logical vector of indices
	des <- cbind(1,out_test)
		
	lmOb <- limma::lmFit(dat_test,design=des)
	lmOb2 <- limma::eBayes(lmOb)
	stats <- lmOb2$t[,2]
	
	cov_mat <- cor(t(dat_test[which(inds1),]))
	
#	if(length(stats)!=(dim(gsea_mat)[1])) print("STOPSTOPSTOPSTOP")
	if(length(stats)!=length(inds1)) print("STOPSTOPSTOPSTOP")
	
	##if(solve(cov_mat)){}
	if(kappa(cov_mat)>1500) {
		increa <- svd(cov_mat)$d[1]/1500 # 1500 is the goal condition number
		cov_mat <- cov_mat + diag(increa,dim(cov_mat)[1])
		}	
	
	test_stat_cor <- t(stats[which(inds1)])%*%solve(cov_mat)%*%(stats[which(inds1)])
	
	te_null_cor <- easier_null(t(dat_test),len_inds=length(which(inds1)),its=200,num_nest=1000,out_test)
	
	#easy_cor <- calc_emp_p(te_null_cor,test_stat_cor[1])
	comp_test_p <- calc_emp_p(te_null_cor,test_stat_cor[1])
	
	# SELF-CONTAINED TEST
	self_test_p <- pchisq(test_stat_cor,df = length(which(inds1)),lower.tail=F)	
	list(comp_test_p,self_test_p)
}





go2eg <- as.list(org.Hs.egGO2EG)
#head(go2eg)

go2alleg <- as.list(org.Hs.egGO2ALLEGS)
#head(go2alleg)

eg2go <- as.list(org.Hs.egGO)

inflamm <- go2eg['GO:0006954'][[1]] # inflamm 
immune <- go2eg['GO:0006955'][[1]] # immune
placen <- go2eg['GO:0001890'][[1]] # placenta dev
apop <- go2eg['GO:0006915'][[1]] # apoptosis




# oncogenic, KEGG
p53_path <- keggGet("path:hsa04115") # p53
mapk_path <- keggGet("path:hsa04010") #mapk
tgf_beta_path <- keggGet("path:hsa04350") #tgf_beta
focal_path <- keggGet("path:hsa04510") #focal adhesion 
cell_prog_path <- keggGet("path:hsa04110") # cell progression


# random
chemokine_path <- keggGet("path:hsa04062") # chemokine pathway (not necessarily one of the cancer pathways, not sure why using it


p53_entrez <- kegg_entrez(p53_path)
mapk_entrez <- kegg_entrez(mapk_path)
tgf_beta_entrez <- kegg_entrez(tgf_beta_path)
focal_entrez <- kegg_entrez(focal_path)
cell_prog_entrez <- kegg_entrez(cell_prog_path)
chemokine_entrez <- kegg_entrez(chemokine_path)


x_gene <- org.Hs.egGENENAME
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- mappedkeys(x_gene)
# Convert to a list
xx <- as.list(x_gene[mapped_genes])

eg_sym <- as.list(org.Hs.egSYMBOL)
sym_eg <- as.list(org.Hs.egALIAS2EG)





geo <- getGEO('GDS3929')
#class(geo@gpl@data.table)

data <- data.table(dataTable(geo)@table)
## data <- dataTable(geo)@table

out <- dataTable(geo)@columns
typ_blood <- as.numeric(dataTable(geo)@columns$specimen)
thresh <- c(7233,8036,8156)




tmp <- sym_eg[(as.character(data[[2]]))]
x <- sapply(tmp,function(x){as.character(x[1])})
na_inds <- is.na(names(x))
#	na_inds <- unlist(sapply(x,function(x){length(x)==0}))
vals <- c(unlist(x),rep(-5,sum(na_inds))) # -5 since won't map to an entrez id 
inds <- c(which(!na_inds),which(na_inds))
dat2_vals_fin <- vals[order(inds)]	


# # ## methods to use:
# # PAGE (citation?), CAMERA (citation?), ROAST, GSEA, MEACA?, 'gsea made simple' by irizarry
# # SAFE? (citation?)
# # GAGE? (citation?)
# # MAST






smoke <- as.numeric((dataTable(geo)@columns)$disease.state) - 1

#data2 <- data[,-c(1,2),with=F]
data2_pre <- data[,-c(1,2),with=F]

design <- model.matrix(~ smoke)
fit <- lmFit(data2_pre, design=design)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, genelist=as.character(unlist(data[,2,with=F])))
tt


thresh <- c(7233,8036,8156)

gsea_list <- lapply(list('p53'=p53_entrez,'mapk'=mapk_entrez,'tgf_beta'=tgf_beta_entrez,'cell_prog'=cell_prog_entrez,'plac'=placen,'inflamm'=inflamm,'immune'=immune),inds_gdc)
gsea_mat <- as.matrix(as.data.frame(gsea_list))


sav_res <- list()

for(type in 1:length(unique(typ_blood))){
#for (type in 1:3){
#type <- 1
	#data2_pre <- data[,-c(1,2),with=F][,which(typ_blood==type),with=F]
	
	row_me <- order(rowMeans(data2_pre[,which(typ_blood==type),with=F]))
	inds_detectable <- tail(row_me,thresh[type])
	data2 <- data2_pre[inds_detectable]
	
	gage_mat <- as.matrix(data2)
	rownames(gage_mat) <- as.character(data[[2]][inds_detectable])
	
# # #(24526 - thresh[type]) : 24
	# row_me <- order(rowMeans(data2_pre))
	# data2 <- data2_pre[,which(typ_blood==type),with=F]
	
	sav_res[[type]] <- list()

	safe_smoke <- safe::safe(as.matrix(data2)[,which(typ_blood==type)],smoke[which(typ_blood==type)],C.mat=gsea_mat[inds_detectable,],Pi.mat=1000)
	
	for (gen_ind in 1:length(gsea_list)){
	
#kk <- 2 
#type <- kk
	camera_rank_smoke <- limma::camera(data2[,which(typ_blood==type),with=F],index=list(gsea_list[[gen_ind]][inds_detectable]) ,design=cbind(1,smoke[which(typ_blood==type)]),use.ranks=T,inter.gene.cor=NA,sort=FALSE)
	camera_norank_smoke <- limma::camera(data2[,which(typ_blood==type),with=F],index=list(gsea_list[[gen_ind]][inds_detectable]),design=cbind(1,smoke[which(typ_blood==type)]),use.ranks=F,inter.gene.cor=NA,sort=FALSE)
	
	#limma::cameraPR(data2, inds_gdc(inflamm),design=smoke)
	roast_smoke <- limma::roast(y= data2[,which(typ_blood==type),with=F], index=list(gsea_list[[gen_ind]][inds_detectable]),design=cbind(1,smoke[which(typ_blood==type)]),sort='none')
	
	gage_test <- gage::gage(gage_mat[,which(typ_blood==type)], gsets = list(first = as.character(data[[2]][inds_detectable][which(gsea_list[[gen_ind]][inds_detectable])])),ref = which(smoke[which(typ_blood==type)]==0),same.dir=F,samp = which(smoke[which(typ_blood==type)]==1),use.fold=F,compare = 'unpaired')
	
	gage_test_rank <- gage::gage(gage_mat[,which(typ_blood==type)], gsets = list(first = as.character(data[[2]][inds_detectable][which(gsea_list[[gen_ind]][inds_detectable])])),ref = which(smoke[which(typ_blood==type)]==0),same.dir=F,samp = which(smoke[which(typ_blood==type)]==1),use.fold=F,rank.test=T,compare = 'unpaired')
	
	# gage_mat[which(gsea_list[[gen_ind]][inds_detectable]),which(typ_blood==type)]
	
	poly_ret <- polymorph_test(gage_mat[,which(typ_blood==type)],smoke[which(typ_blood==type)],gsea_mat[inds_detectable,gen_ind])
	
	comp_tes <- poly_ret[[1]]
	self_tes <- poly_ret[[2]]
	
	tm_lis <- list(camera_rank_smoke[[4]][1],
	camera_norank_smoke[[4]][1],
	roast_smoke[[5]][1],
	safe_smoke@global.pval[gen_ind],   ## safe seems to already be somewhat rank-based
	gage_test$greater[3],
	gage_test_rank$greater[3],
	comp_tes,
	self_tes
	)
	sav_res[[type]][[gen_ind]] <- tm_lis
	print(gen_ind)
	
	}
		
	print(type)
}


#save(sav_res,file='~/Google_Drive/Oslo_sync/gsea_paper/gsea_analysis_ob_detectable.RObject',compress=F)
#save(sav_res_sav,file='~/Google_Drive/Oslo_sync/gsea_paper/gsea_analysis_ob_detectable_GAGEfix.RObject',compress=F)



# ## FOR FIXING THE ISSUE WITH GAGE
# for (type in 1:3){
	# for (gen_ind in 1:7){
		# tmp <- sav_res_sav[[type]][[gen_ind]]
		
		# tmp[[5]] <- sav_res[[type]][[gen_ind]][[5]]
		# tmp[[6]] <- sav_res[[type]][[gen_ind]][[6]]
		
		# sav_res_sav[[type]][[gen_ind]][[5]] <- tmp[[5]]
		# sav_res_sav[[type]][[gen_ind]][[6]] <- tmp[[6]]
		# #sav_res[[type]][[gen_ind]] <- tm_lis
		# }	
	# }


col_nams <- names(gsea_list)
row_nams <- c('camera_rank_smoke','camera_norank_smoke','roast_smoke','safe_smoke','gage_test','gage_test_rank','competitive_test','self_cont_test')

unlis_sav_res <- unlist(sav_res)
# unlis_sav_res <- unlist(sav_res_sav2)
# unlis_sav_res <- unlist(sav_res_sav)

len3 <- length(unlis_sav_res)/3

# rep(c(T,F,F),each=len3)
# rep(c(F,T,F),each=len3)
# rep(c(F,F,T),each=len3)

# mat_per2 <- matrix(unlis_sav_res[rep(c(T,F,F),each=len3)],nrow=length(row_nams)) # type maternal peripheral
# neonat_cord2 <- matrix(unlis_sav_res[rep(c(F,T,F),each=len3)],nrow=length(row_nams)) # type neonatal cord
# plac2 <- matrix(unlis_sav_res[rep(c(F,F,T),each=len3)],nrow=length(row_nams)) # type placenta

mat_per <- matrix(unlis_sav_res[rep(c(T,F,F),each=len3)],nrow=length(row_nams)) # type maternal peripheral
neonat_cord <- matrix(unlis_sav_res[rep(c(F,T,F),each=len3)],nrow=length(row_nams)) # type neonatal cord
plac <- matrix(unlis_sav_res[rep(c(F,F,T),each=len3)],nrow=length(row_nams)) # type placenta

# identical(plac,plac2)
# identical(neonat_cord,neonat_cord2)
# identical(mat_per,mat_per2)

colnames(mat_per) <- col_nams
colnames(neonat_cord) <- col_nams
colnames(plac) <- col_nams

rownames(mat_per) <- row_nams
rownames(neonat_cord) <- row_nams
rownames(plac) <- row_nams

library(xtable)

xtable(mat_per)
xtable(neonat_cord)
xtable(plac)






### TRYING TO CONFIRM RESULTS OF THE VOTAVOVA PLACENTA PAPER FROM WHICH WE GET THE DATA, AND INDEED SEEM TO

# periph: 7233
# plac 8156
# cord : 8036 

typ_blood <- as.numeric(dataTable(geo)@columns$specimen)

# 1 = maternal peripheral
# 2 = neonatal cord
# 3 = term placenta

# thresholds for detection
thresh <- c(7233,8036,8156)


for ( type in 1:3){
#type <- 1
data2_pre <- data[,-c(1,2),with=F][,which(typ_blood==type),with=F]

#(24526 - thresh[type]) : 24
row_me <- order(rowMeans(data2_pre))
inds_detectable <- tail(row_me,thresh[type])
data2 <- data2_pre[inds_detectable]

design <- model.matrix(~ smoke[which(typ_blood==type)])
fit <- lmFit(data2, design=design)
fit <- eBayes(fit)
#tt <- topTable(fit, coef=2,number=20, genelist=as.character(unlist(data[,2,with=F])))
tt <- topTable(fit, coef=2,number=20, genelist=as.character(unlist(data[inds_detectable,2,with=F])))

print(tt[,1])
print("")
print("")
print("")
print("")
}



