#' Form graph
#'
#' @param path filepath prefix
#' @param edge_func function to determine if there's an edge
#' @param K number of nodes or vector of indicies
#' @param max_cluster maximum number of clusters
#' @param cores number of cores for parallelization
#' @param verbose boolean
#'
#' @return \code{igraph} object
#' @export
form_graph <- function(path = "/raid6/Kevin/supernodeWGS_results/blocks/",
                       edge_func = function(x){abs(mean(as.numeric(unlist(x))))>0.1},
                       K = 200, max_cluster = NA, cores = 1, verbose = T){
  doMC::registerDoMC(cores = cores)
  if(length(K) == 1){
    max_cluster <- K
    combn_mat <- utils::combn(K, 2)
    idx_mat <- combn_mat
  } else {
    stopifnot(!is.na(max_cluster))
    num_node <- length(K)
    idx_mat <- utils::combn(num_node, 2)
    combn_mat <- rbind(K[idx_mat[1,]], K[idx_mat[2,]])
  }

  func <- function(i){
    if(verbose & i %% floor(ncol(combn_mat)/10) == 0) print('*')
    idx <- pair_to_index(combn_mat[1,i], combn_mat[2,i], max_val = max_cluster)

    cor_block <- 0 #debugging reasons
    load(paste0(path, idx, ".RData"))
    edge_func(cor_block)
  }

  i <- 0 #debugging reasons
  res <- foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), func(i))
  idx <- which(as.numeric(res) == 1)
  edges <- idx_mat[,idx]

  g <- igraph::graph.empty(n = num_node, directed = F)
  if(length(K) > 1){
    igraph::V(g)$name <- as.character(K)
  }

  igraph::add_edges(g, edges = edges)
}

#' Form correlation
#'
#' @param path filepath prefix
#' @param summary_func function to summarize each edge
#' @param K number of nodes or vector of indicies
#' @param max_cluster maximum number of clusters
#' @param cores number of cores for parallelization
#' @param verbose boolean
#'
#' @return list
#' @export
form_correlation <- function(path = "/raid6/Kevin/supernodeWGS_results/blocks/",
                             summary_func = function(x){mean(as.numeric(unlist(x)))},
                             K = 200, max_cluster = NA, cores = 1, verbose = T){
  doMC::registerDoMC(cores = cores)
  if(length(K) == 1){
    max_cluster <- K
    num_node <- K
    combn_mat <- utils::combn(num_node, 2)
  } else {
    stopifnot(!is.na(max_cluster))
    idx_mat <- utils::combn(length(K), 2)
    combn_mat <- rbind(K[idx_mat[1,]], K[idx_mat[2,]])
  }

  func <- function(i){
    if(verbose & i %% floor(ncol(combn_mat)/10) == 0) print('*')
    idx <- pair_to_index(combn_mat[1,i], combn_mat[2,i], max_val = max_cluster)

    cor_block <- 0 #debugging reasons
    load(paste0(path, idx, ".RData"))
    summary_func(cor_block)
  }

  i <- 0 #debugging reasons
  foreach::"%dopar%"(foreach::foreach(i = 1:ncol(combn_mat)), func(i))
}

#' Form test statistics
#'
#' @param vec vector of test statistics
#' @param clustering vector of clustering assignments
#' @param flag_vec boolean vector of same length as \code{vec}, where indices with
#' \code{TRUE} will be omitted from calculation
#' @param path filepath
#' @param K number of nodes or vector of indicies
#' @param max_cluster maximum number of clusters
#' @param sparse boolean for sparse PCA
#' @param sumabsv paramter for sparse PCA
#' @param cores number of cores for parallelization
#' @param verbose boolean
#' @param ... extra parameters for PCA::SPC
#'
#' @return vector
#' @export
form_testvec <- function(vec, clustering, flag_vec = rep(FALSE, length(vec)),
                      path = "/raid6/Kevin/supernodeWGS_results/blocks/",
                      K = 200, max_cluster = NA, sparse = F, sumabsv = 4, cores = 1, verbose = T){
  stopifnot(length(vec) == length(clustering))
  doMC::registerDoMC(cores = cores)

  if(length(K) == 1) {
    cluster_vec <- 1:K
    max_cluster <- K
  } else {
    cluster_vec <- K
    stopifnot(!is.na(max_cluster))
  }

  func <- function(i){
    if(verbose & i %% floor(max(K)/10) == 0) print('*')

    cluster_idx <- which(clustering == i)
    testvec <- vec[cluster_idx]; testflag <- flag_vec[cluster_idx]
    keep_idx <- intersect(which(!is.na(testvec)), which(!testflag))
    testvec <- testvec[keep_idx]

    idx <- pair_to_index(i, i, max_val = max_cluster)
    cor_block <- 0 #debugging reasons
    load(paste0(path, idx, ".RData"))
    cor_block <- as.matrix(cor_block)
    cor_block <- cor_block[keep_idx, keep_idx, drop = F]
    if(length(cor_block) == 0) return(list(val = NA, index = NA, sparsity = NA))

    eig <- .extract_eigenvector(cor_block, sparse = sparse, sumabsv = sumabsv)

    unsigned <- sqrt(as.numeric((eig%*%testvec)^2 / t(eig) %*% cor_block %*% eig))
    sign <- .determine_sign(eig %*% t(eig) %*% testvec)

    stopifnot(length(eig) == length(keep_idx))

    list(val = sign * unsigned, index = cluster_idx[keep_idx[which(abs(eig) >= 1e-5)]],
         sparsity = length(which(abs(eig) >= 1e-5))/length(eig))
  }

  i <- 0 #debugging reasons
  res <- foreach::"%dopar%"(foreach::foreach(i = cluster_vec), func(i))

  vec <- sapply(res, function(x){x$val})
  index <- lapply(res, function(x){x$index})
  sparsity <- sapply(res, function(x){x$sparsity})

  list(vec = vec, index = index, sparsity = sparsity)
}

.determine_sign <- function(vec){
  pos <- length(which(vec > 0))
  neg <- length(which(vec < 0))
  if(pos > neg) return(1)
  if(neg > pos) return(-1)

  stats::rbinom(1, 1, 0.5)*2-1
}

.extract_eigenvector <- function(mat, sparse = F, sumabsv = 4){
  if(sparse & ncol(mat) > sumabsv^2){
    eig <- PMA::SPC(mat, trace = F, sumabsv = sumabsv)
    as.numeric(eig$v[,1])
  } else {
    eig <- eigen(mat)
    as.numeric(eig$vectors[,1])
  }
}

