##########################################################################
## An et al. 2018
## De novo risk score analysis
##
## The analysis is consist of the following steps:
## 1. Limit the analysis to the set of rare categories 
## 2. Fit a Lasso model using the 519 families 
##    previously described in Werling et al. 2018,
##    and compute the predictive R2 on the remaining 1,383 new families
##########################################################################

rm(list=ls())
gc()
library(data.table)
library(glmnet)


###############################################
## 0. Load data
###############################################

## Deonvo data
adj.denovo.counts <- fread(
	"gunzip -c result.sumVar.20171115_allSets_allAF_cwasSim_20180108.trimmed.adjusted_denovo_score.txt.gz",
	data.table=FALSE)
row.names(adj.denovo.counts) <- adj.denovo.counts[, 1]
adj.denovo.counts[, 1] <- NULL
  
## We may consider a subset of categories,
## for example, "nonredundant", "noncoding", "promoter", "coding".
## This .txt file contains the categories to be included in this analysis.
list.subset = fread("list.promoter.txt", header=TRUE, data.table=FALSE)
i.match <- na.omit(match(list.subset[, 1], colnames(adj.denovo.counts)))
adj.denovo.counts <- adj.denovo.counts[, i.match]

  
#################################################
## 1. Limit the analysis to the set of rare categories 
##
## In particular, we exclude de novo events falling in an annotation category 
## if their total frequency, adjusted for paternal age, 
## taken over an annotation category, 
## is >= 3 in controls.
#################################################

## data with de novo counts adjusted for paternal age
nrd.1902 <- fread("gunzip -c result.perm_p.20171115_allSets_allAF_cwasSim_20180108.2.txt.gz",
                  sep="\t", data.table=FALSE)

## limit to rare events
ct.thres=3
annot.rare <- nrd.1902$Annotation_combo[which(nrd.1902$Con_count_adj < ct.thres)]
i.screen = na.omit(match(annot.rare, colnames(adj.denovo.counts)))
covariates = as.matrix(adj.denovo.counts[, i.screen])

#################################################
## 2. Fit a Lasso model
##
## The 519 families are used as training set, 
## and the 1,383 new families as testing set.
## Response are -1 for controls and 1 for cases.
#################################################

## Get family ID and sample ID of all 1902 families
fam.1902 <- data.frame(strsplit(row.names(covariates), "_"))
fam.1902 <- data.frame(t(fam.1902))
colnames(fam.1902) <- c("familyID", "sampleID")

## the list of previous 519 families
fam.519 <- fread("infoForDeNovoRateCorrection-170320.txt", header=TRUE)
i.519 <- which(fam.1902$familyID %in% fam.519$familyID)

## the 1,383 new families as testing set
test.covariates = covariates[-i.519, ]
i.case.test <- grep("_p", row.names(test.covariates))
i.control.test <- grep("_s", row.names(test.covariates))
response.test <- rep(-1, nrow(test.covariates))
response.test[i.case.test] <- 1

## the 519 families as training set
covariates = covariates[i.519, ]
N <- nrow(covariates)
i.case <- grep("_p", row.names(covariates))
i.control <- grep("_s", row.names(covariates))
response <- rep(-1, N)
response[i.case] <- 1

## fit lasso 
## use cross validation to choose optimal lambda
fold = 5
n.family = N/2
n.test = round(n.family / fold)
seeds = c(99, 109, 119, 129, 139, 149, 159, 169, 179, 189)
results = vector("list", length(seeds))

for (i in c(1:length(seeds))) {
  set.seed(seeds[i])

  ## randomly split data into 5 folds for cross validation
  i.permute = sample(c(1:n.family), size=n.family, replace = FALSE)
  foldid = rep(0, N)
  for (fo in 1:fold){
    if (fo < fold) {
      i.family.test = i.permute[((fo-1)*n.test + 1) : (fo*n.test)]
    } else {
      i.family.test = i.permute[((fo-1)*n.test + 1) : n.family]
    }
    i.test = sort(c(2*i.family.test, 2*i.family.test-1))
    foldid[i.test] = fo
  }
  
  ## cross-validation to choose lamda
  cv.out = cv.glmnet(covariates, 
                     response, 
                     foldid=foldid, 
                     alpha=1, type.measure="deviance",
                     standardize=TRUE, nlambda=20)
  
  lambda.choose = cv.out$lambda.min
  i.choose = match(lambda.choose, cv.out$lambda)
  
  ## get lasso coefficient
  full.glmnet = cv.out$glmnet.fit
  lasso.betas = as.matrix(full.glmnet$beta)[, i.choose]
  n.select = sum(abs(lasso.betas) > 0)
  
  ## predictive R2 on held out set
  predict.y = predict(full.glmnet, newx=test.covariates)[, i.choose]
  mse = mean((predict.y - response.test)^2)
  R_sq = 1 - mse

  ## save
  results[[i]] = list(R_sq = R_sq,
  					  lambda.choose=lambda.choose, 
                      lasso.betas=lasso.betas,                 
                      i.choose=i.choose, 
                      n.select=n.select)
}

## save results
save(results, seeds, file="out_predict_r2.RData")

## average predictive R2
R_sqs = sapply(results, function(x){x$R_sq})
print(paste("Average predictive R2 is", mean(R_sqs)))



