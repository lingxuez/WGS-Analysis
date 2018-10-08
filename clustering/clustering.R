#########################################################
## An et al. 2018
## Spectral clustering of categories
##
## This analysis includes the following steps:
## 1. Compute the correlation structure
## 2. Estimate the top eigenvectors
## 3. Use t-sne projection to initialize cluster centers
## 4. Apply spectral clustering
#########################################################

rm(list=ls())
gc()

library(data.table)
library(doParallel)
library(matrixStats)
library(HiClimR)
library(rARPACK)
library(tsne)
library(fastcluster)
library(mclust)
library(cluster)


#########################################################
## 1. Compute correlation structure
##
## We transformed p-values to z-scores.
## Because of zero counts and p-values of one, 
## we added one pseudo-count to both case and control counts.
#########################################################
## Load the count matrix for probands and siblings;
## each row is one simulation, and each column is one annotation
pro_counts = fread("zcat table.Pro_count.fam1902_20180117.txt.gz",
                   drop=1, sep="\t", header=TRUE, data.table=FALSE)
sib_counts = fread("zcat table.Sib_count.fam1902_20180117.txt.gz"),
                   drop=1, sep="\t", header=TRUE, data.table=FALSE)

## limit to promoters 
list.subset = fread("list.promoter.txt",
                    header=TRUE, data.table=FALSE)
i.keep = match(list.subset[, 1], colnames(pro_counts))
print(paste(length(i.keep), "annotations are retained for analysis."))

## one-sided z-value after adding one pseudo count
get_zscore = function(x, y) {
  ## x: counts in probands; y: counts in siblings
  p_values = sapply(1:length(x), function(i){
    binom.test(x[i]+1, x[i]+y[i]+2, p=0.5, "greater")$p.value
    })
  ## convert to z-scores
  return(-qnorm(as.numeric(p_values)))
}

## compute z-scores for all annotations in parallel
mc.cores = 10
zscores = mcmapply(function(j){get_zscore(pro_counts[, j], sib_counts[, j])},
               i.keep, mc.cores=mc.cores)
colnames(zscores) = colnames(sib_counts)[i.keep]


#########################################################
## 2. Estimate the top eigenvectors
#########################################################
nSims = nrow(zscores)
nAnnot = ncol(zscores)

## compute the correlation matrix
corr.pseudo = HiClimR::fastCor(zscores, nSplit=2)

## compute the negative Laplacian matrix to normalize node degrees
neg_Laplacian = abs(corr.pseudo)
degrees = colSums(neg_Laplacian)
for (j in 1:nrow(corr.pseudo)) {
  neg_Laplacian[j, ] <- neg_Laplacian[j, ] / sqrt(degrees)
  neg_Laplacian[, j] <- neg_Laplacian[, j] / sqrt(degrees)
}

## It suffices to get the first Nmin eigenvalues and eigenvectors
## because the rank of the negative Laplacian is at most Nmin
Nmin = min(nSims, nAnnot)

lap.eigs = rARPACK::eigs(neg_Laplacian, Nmin)
eigs.values = lap.eigs$values
eigs.vectors = lap.eigs$vectors

## save the results
fwrite(as.data.frame(eigs.values), 
       file="negLaplacian_eigValues_Z_1side_pseudo.txt", 
       sep="\t", col.names=FALSE)
fwrite(as.data.frame(eigs.vectors), 
       file=p"negLaplacian_eigVectors_Z_1side_pseudo.txt", 
       sep="\t", col.names=FALSE)

#########################################################
## 3. Use t-sne projection to initialize cluster centers
#########################################################
## How many eigen vectors to use in spectral clustering
num.dim = 50

## ignored the first eigenvector because it mainly accounts for the mean level
U.pick = eigs.vectors[, 2:(num.dim)]
## normalized the remaining 49 eigenvectors
U.norm = t(scale(t(U.pick), center=FALSE, scale=sqrt(rowSums(U.pick^2))))

## tsne
tsne.out = tsne(U.norm, k=2, perplexity=30, max_iter=500, whiten=TRUE)
tsne.out = as.data.frame(tsne.out)
colnames(tsne.out) = c("tsne1", "tsne2")

#########################################################
## 4. Apply spectral clustering
#########################################################
K = 70

## fit K-means on tsne for center initialization
km_tsne = kmeans(tsne.out, centers=K, nstart=300, iter.max=100)
  
## find the closest point to each center and use to initialize K-means
dist.to.centor.tsne = sqrt(rowSums((as.matrix(tsne.out) - 
                                    km_tsne$centers[km_tsne$cluster, ])^2))
i.init.pt = sapply(1:K,
                   function(k){
                    i.check = which(km_tsne$cluster == k)
                    i.min = which.min(dist.to.centor.tsne[i.check])
                    return(i.check[i.min])
                    })
  
## K-means
## no re-start because initial centers are given
fit = kmeans(U.norm, centers=U.norm[i.init.pt, ], iter.max=iter.max)
  
## results
save(fit, U.norm, file="out_clustering.RData")




