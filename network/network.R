## The following script is for the network analysis, which follows the clustering analysis.

rm(list=ls())
library(supernodeWGS) #make sure to install the supernodeWGS R package first, provided.

suffix <- "promoter"
count_threshold <- 20
cor_threshold <- 0.12
size_threshold <- 2
cores <- 4

#some cleanup first
### load in the appropriate clustering
load(paste0("supernodeWGS_results/clustering_", suffix, ".RData"))
dat1 <- clustering
K <- max(dat1$Cluster)
clustering <- dat1$Cluster
name1 <- as.character(dat1$Annotation)

dat2 <- read.csv(paste0("supernodeWGS/data/result.perm_p.20171115_allSets_allAF_cwasSim_20180108.", suffix, ".txt"), sep = "\t")
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

cor_mat <- supernodeWGS::form_correlation(path = paste0("supernodeWGS_results/blocks_", suffix, "/"),
                            K = cluster_idx, max_cluster = K, cores = cores)

g <- supernodeWGS::form_graph_from_correlation(cor_mat, func = function(x){x>cor_threshold}, K = cluster_idx)
adj <- as.matrix(igraph::as_adjacency_matrix(g))

########
testvec_res <- supernodeWGS::form_testvec(vec1, clustering, flag_vec,
                                  path = paste0("supernodeWGS_results/blocks_", suffix, "/"),
                                  K = cluster_idx, max_cluster = K, cores = cores, sparse = T,
                                  sumabsv = 5.25)
zval_supernode <- testvec_res$vec
annotation_cluster <- testvec_res$index

testvec_res2 <- supernodeWGS::form_testvec(vec1, clustering, flag_vec,
                                             path = paste0("supernodeWGS_results/blocks_", suffix, "/"),
                                             K = cluster_idx, max_cluster = K, cores = cores)

riskvec_res <- supernodeWGS::form_testvec(risk1, clustering, flag_vec2,
                                             path = paste0("supernodeWGS_results/blocks_", suffix, "/"),
                                             K = cluster_idx, max_cluster = K, cores = cores)
tmp <- riskvec_res$vec
tmp[is.na(tmp)] <- 0
risk_supernode <- exp(tmp)

##########

res <- supernodeWGS::hmrf(zval_supernode, adj,
                          seedindex = rep(0, length(cluster_idx)), verbose = T)
fdr <- supernodeWGS::report_results(cluster_idx, 1-res$post, 1-stats::pnorm(zval_supernode), res$Iupdate)

##########

# dawn graph, with color based on the significance
node_size <- sapply(cluster_idx, function(x){
  length(which(clustering == x))
})
node_col <- rep("white", length(zval_supernode))
set.seed(10)
igraph::plot.igraph(g, vertex.size = node_size, vertex.label = NA,
                    vertex.color = node_col, main = "Gene graph")

#######

### output CSVs
# pvalues
zz <- rep(NA, length(zval_supernode))
zz[which(risk_supernode >= 1)] <- 1-pnorm(zval_supernode[which(risk_supernode >= 1)])
zz[which(risk_supernode < 1)] <- pnorm(zval_supernode[which(risk_supernode < 1)])

mat <- data.frame(Cluster.idx = as.numeric(igraph::V(g)$name), Pval.cluster = zz, Risk = risk_supernode)
write.csv(mat, file = paste0("../results/pvalue_risk",
                             suffix, ".csv"), quote = F)

# graph

adj2 <- as.matrix(igraph::as_adj(g))

stopifnot(all(cluster_idx == as.numeric(igraph::V(g)$name)))

write.table(adj, file = paste0("../results/graph_", suffix, ".csv"), sep = ",", row.names = T)

# graph layout
set.seed(10)
maxz <- max(zval_supernode)
minz <- min(zval_supernode)

node_size <- sapply(cluster_idx, function(x){
  length(which(clustering == x))
})
node_size <- (node_size - min(node_size))/max(node_size)
node_size <- 2*(node_size+1)^2

node_col <- sapply(zval_supernode, function(x){
  tmp <- (maxz-x)/(maxz-minz)
  tmp <- 1/(1+exp(-10*(tmp - 0.2)))
  rgb(1, tmp, tmp)
})
node_col[1] <- rgb(0,1,0)
node_col[7] <- rgb(0,0,1)

l <- igraph::layout_nicely(g)
igraph::plot.igraph(g, layout = l, vertex.size = node_size, vertex.label = NA,
                    vertex.color = node_col, main = "Pre-DAWN")

layout_mat <- as.data.frame(cbind(as.numeric(igraph::V(g)$name), l))
names(layout_mat) <- c("Cluster.index", "X.pos", "Y.pos")
write.csv(layout_mat, file = paste0("../results/graph_layout_",
                                    suffix, ".csv"), quote = F)

# clusterings
mode = "all" #(mode can be "all" or "fdr")

load(paste0("../results/results_", suffix, ".RData"))

if(mode == "fdr"){
  tmp_idx <- which(fdr$FDR <= 0.0105)
} else {
  tmp_idx <- 1:igraph::vcount(g)
}

cluster_indices <- as.numeric(igraph::V(g)$name)[tmp_idx]
stopifnot(all(sort(cluster_indices) == sort(cluster_idx[tmp_idx])))
results <- vector("list", length(cluster_indices))

for(i in 1:length(cluster_indices)){
  node_idx <- testvec_res$index[[tmp_idx[i]]]
  node_nam <- dat1$Annotation[node_idx]

  stopifnot(length(unique(clustering[node_idx])) == 1)
  stopifnot(length(node_idx) <= size_vec[cluster_indices[i]])

  results[[i]] <- node_nam
}
len <- max(sapply(results, length))

mat <- matrix(" ", nrow = len, ncol = length(cluster_indices))
for(i in 1:length(cluster_indices)){
  mat[1:length(results[[i]]),i] <- as.character(results[[i]])
}
colnames(mat) <- cluster_indices
mat <- as.data.frame(mat)

write.csv(mat, file = paste0("../results/cluster_annotations_tableformat_",
                             suffix, "_", mode, ".csv"), quote = F)



