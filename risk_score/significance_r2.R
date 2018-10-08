##########################################################################
## An et al. 2018
## P-value of predictive R2
##
## We use permutation to estimate the p-value of predictive R2.
## In each permutation, we randomly flipped the labels of case and control
## with probability 0.5 in each family.
## Next, we repeated the full procedure, including:
## 1. Limit the analysis to the set of rare categories 
## 2. Fit a Lasso model using the previously described 519 families 
##    and compute the predictive R2 on the remaining 1,383 new families
##########################################################################

rm(list=ls())
gc()
library(data.table)
library(glmnet)


###############################################
## Load Data
###############################################

## Deonvo data
adj.denovo.counts <- fread(
  "gunzip -c result.sumVar.20171115_allSets_allAF_cwasSim_20180108.trimmed.adjusted_denovo_score.txt.gz",
  data.table=FALSE)
row.names(adj.denovo.counts) <- adj.denovo.counts[, 1]
adj.denovo.counts[, 1] <- NULL

## Get family ID and sample ID of all 1902 families
fam.1902 <- data.frame(strsplit(row.names(adj.denovo.counts), "_"))
fam.1902 <- data.frame(t(fam.1902))
colnames(fam.1902) <- c("familyID", "sampleID")

## We may consider a subset of categories,
## for example, "nonredundant", "noncoding", "promoter", "coding".
## This .txt file contains the categories to be included in this analysis.
list.subset = fread("list.promoter.txt", header=TRUE, data.table=FALSE)
i.match <- na.omit(match(list.subset[, 1], colnames(adj.denovo.counts)))
adj.denovo.counts <- adj.denovo.counts[, i.match]

## the list of previous 519 families
fam.519 <- fread("infoForDeNovoRateCorrection-170320.txt", header=TRUE)
i.519 <- which(fam.1902$familyID %in% fam.519$familyID)


##########################
## Permutation Test
###########################
npermute = 1000
ct.thres = 3

original.id = row.names(adj.denovo.counts)
results = vector("list", npermute)

for (i in c(1:npermute)) {
  
  ######################
  ## 0. Randomly flipped the labels of case and control in each family
  ######################
  n.fam = length(original.id)/2
  permute.id = original.id
  for (fa in c(1:n.fam)) {
    if (runif(1) < 0.5) {
      ## flip
      permute.id[2*fa-1] = original.id[2*fa]
      permute.id[2*fa] = original.id[2*fa-1]
    }
  }
  
  ######################
  ## 1. Limit the analysis to the set of rare categories 
  ######################
  i.pro.total <- grep("_p", permute.id)
  i.sib.total <- grep("_s", permute.id)

  ## rare events
  total_count = colSums(adj.denovo.counts[i.pro.total, ])
  annot.rare <- colnames(adj.denovo.counts)[which(total_count < ct.thres)]
  i.screen = na.omit(match(annot.rare, colnames(adj.denovo.counts)))
  covariates = as.matrix(adj.denovo.counts[, i.screen])
  
  #######################
  ## 2. Train Lasso model using 519 families
  #######################
  test.covariates = covariates[-i.519, ]
  i.case.test <- grep("_p", permute.id[-i.519])
  i.control.test <- grep("_s", permute.id[-i.519])
  
  response.test <- rep(-1, nrow(test.covariates))
  response.test[i.case.test] <- 1
  
  covariates = covariates[i.519, ]
  N <- nrow(covariates)
  i.case <- grep("_p", permute.id[i.519])
  i.control <- grep("_s", permute.id[i.519])
  
  response <- rep(-1, N)
  response[i.case] <- 1
  
  ## use cross validation to choose lambda for lasso
  ## randomly split 519 families to 5 folds
  fold = 5
  n.family = N/2
  n.test = round(n.family / fold)
  
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
                     standardize=TRUE)
  
  ## final lasso model
  lambda.choose = cv.out$lambda.min
  i.choose = match(lambda.choose, cv.out$lambda)
  full.glmnet = cv.out$glmnet.fit
  lasso.betas = as.matrix(full.glmnet$beta)[, i.choose]
  n.select = sum(abs(lasso.betas) > 0)
  
  ## R sq on the testing set
  predict.y = predict(full.glmnet, newx=test.covariates)[, i.choose]
  mse = mean((predict.y - response.test)^2)
  R_sq = 1 - mse
  
  ## save results
  results[[i]] = list(R_sq = R_sq,
                      lambda.choose=lambda.choose, 
                      i.choose=i.choose, 
                      n.select=n.select
                      )
  
}
  

## save
save(results, file="out_permute_r2.RData"))

## use mean and std to estimuate the null distribution
R_sqs = sapply(results, function(x){x$R_sq})
print(paste("Estimated null distribution: Normal with mean", mean(R_sqs), ", std", sd(R_sqs)))



