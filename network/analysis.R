rm(list=ls())
library(supernodeWGS)

# suffix <- "coding"
suffix <- "promoter"
# suffix <- "noncoding"
count_threshold <- 20
# con_count <- 3
if(suffix == "promoter") {
  cor_threshold <- 0.12 #for promoter
} else if(suffix == "coding") {
  cor_threshold <- 0.07 #for coding
} else {
  cor_threshold <- 0.05 #for noncoding
}
size_threshold <- 2
cores <- 4

#some cleanup first
load(paste0("/raid6/Kevin/supernodeWGS_results/clustering_", suffix, ".RData"))
dat1 <- clustering
K <- max(dat1$Cluster)
clustering <- dat1$Cluster
name1 <- as.character(dat1$Annotation)

if(suffix != "coding"){
  dat2 <- read.csv(paste0("~/supernodeWGS/data/result.perm_p.20171115_allSets_allAF_cwasSim_20180108.", suffix, ".txt"), sep = "\t")
} else {
  dat2 <- read.csv("~/supernodeWGS/data/result.perm_p.20171115_allSets_allAF_cwasSim_20180108.2.txt", sep = "\t")
}
idx2 <- intersect(which(dat2$Adjusted_relative_risk <= 1), which(dat2$Perm_p < 0.5))
vec2 <- dat2$Perm_p
vec2[idx2] <- 1 - vec2[idx2]
risk2 <- dat2$Adjusted_relative_risk

#match up the indices
name2 <- as.character(dat2$Annotation_combo)
idx <- supernodeWGS::matching(name1, name2)
vec1 <- vec2[idx]
risk1 <- risk2[idx]

#########

# find out which indices have too few permutations, and then compute test statistics
count_vec2 <- dat2$Total_count_raw
count_vec1 <- count_vec2[idx]
flag_vec <- rep(FALSE, length(count_vec1))
flag_vec[which(count_vec1 <= count_threshold)] <- TRUE

# #keep only rare variants
# con_adj2 <- dat2$Con_count_adj
# con_adj1 <- con_adj2[idx]
# flag_vec[which(con_adj2 > con_count)] <- TRUE

#add to flag_vec: indicies where vec1 is NA (meaning no match) or screen the names
flag_vec[which(is.na(vec1))] <- TRUE
tab <- supernodeWGS::tabulate_names(name1)
if(suffix != "coding"){
  bool_screen <- supernodeWGS::screen_name(name1)
  flag_vec[which(!bool_screen)] <- TRUE
}

size_vec <- supernodeWGS::cluster_size(vec1, clustering, flag_vec)
cluster_idx <- which(size_vec >= size_threshold)

#convert into z-scores
vec1[!is.na(vec1)] <- stats::qnorm(1 - vec1[!is.na(vec1)])
veclim <- max(abs(vec1[which(!is.infinite(vec1))]), na.rm = T)
vec1[intersect(which(is.infinite(vec1)), which(vec1 > 0))] <- 1.1*veclim
vec1[intersect(which(is.infinite(vec1)), which(vec1 < 0))] <- -1.1*veclim
risk1 <- log(risk1)

flag_vec2 <- flag_vec
flag_vec2[is.infinite(risk1)] <- TRUE

#########

cor_mat <- supernodeWGS::form_correlation(path = paste0("/raid6/Kevin/supernodeWGS_results/blocks_", suffix, "/"),
                            K = cluster_idx, max_cluster = K, cores = cores)

g <- supernodeWGS::form_graph_from_correlation(cor_mat, func = function(x){x>cor_threshold}, K = cluster_idx)
adj <- as.matrix(igraph::as_adjacency_matrix(g))

########
testvec_res <- supernodeWGS::form_testvec(vec1, clustering, flag_vec,
                                  path = paste0("/raid6/Kevin/supernodeWGS_results/blocks_", suffix, "/"),
                                  K = cluster_idx, max_cluster = K, cores = cores, sparse = T,
                                  sumabsv = 5.25)
zval_supernode <- testvec_res$vec
annotation_cluster <- testvec_res$index

testvec_res2 <- supernodeWGS::form_testvec(vec1, clustering, flag_vec,
                                             path = paste0("/raid6/Kevin/supernodeWGS_results/blocks_", suffix, "/"),
                                             K = cluster_idx, max_cluster = K, cores = cores)

riskvec_res <- supernodeWGS::form_testvec(risk1, clustering, flag_vec2,
                                             path = paste0("/raid6/Kevin/supernodeWGS_results/blocks_", suffix, "/"),
                                             K = cluster_idx, max_cluster = K, cores = cores)
tmp <- riskvec_res$vec
tmp[is.na(tmp)] <- 0
risk_supernode <- exp(tmp)

#######

res <- supernodeWGS::hmrf(zval_supernode, adj,
                          seedindex = rep(0, length(cluster_idx)), verbose = T)
fdr <- supernodeWGS::report_results(cluster_idx, 1-res$post, 1-stats::pnorm(zval_supernode), res$Iupdate)

#######

if(suffix == "promoter"){
  threshold_vec <- c(0.07, 0.1, 0.15) # for promoter
} else {
  threshold_vec <- c(0.025, 0.04, 0.06) # for noncoding
}

g_list <- vector("list", length(threshold_vec))
res_list <- vector("list", length(threshold_vec))
fdr_list <- vector("list", length(threshold_vec))

for(i in 1:length(threshold_vec)){
  g_list[[i]] <- supernodeWGS::form_graph_from_correlation(cor_mat, func = function(x){x>threshold_vec[i]}, K = cluster_idx)
  adj_tmp <- as.matrix(igraph::as_adjacency_matrix(g_list[[i]]))

  res_list[[i]] <- supernodeWGS::hmrf(zval_supernode, adj_tmp,
                            seedindex = rep(0, length(cluster_idx)), verbose = T)
  fdr_list[[i]] <- supernodeWGS::report_results(cluster_idx, 1-res_list[[i]]$post, 1-stats::pnorm(zval_supernode), res_list[[i]]$Iupdate)
}

save.image(paste0("/raid6/Kevin/supernodeWGS_results/results_", suffix, ".RData"))
