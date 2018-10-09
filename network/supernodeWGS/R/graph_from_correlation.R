#' Create graph for correlation
#'
#' The input, \code{cor_mat}, follows the format given by \code{form_correlation}.
#' This is important if \code{K} is a vector of indices, as opposed a positive
#' integer for how many nodes there are.
#'
#' @param cor_mat list of correlation values
#' @param func function to determine presence of edge
#' @param K number of nodes or vector of indicies
#'
#' @return graph
#' @export
form_graph_from_correlation <- function(cor_mat, func = function(x){x>0.15}, K = 200){
  if(length(K) == 1){
    stopifnot(length(cor_mat) == K*(K-1)/2)
    idx_mat <- utils::combn(K, 2)
  } else {
    stopifnot(length(cor_mat) == length(K)*(length(K)-1)/2)
    num_node <- length(K)
    idx_mat <- utils::combn(num_node, 2)
  }

  idx <- which(sapply(unlist(cor_mat), func))
  edges <- idx_mat[,idx]

  g <- igraph::graph.empty(n = num_node, directed = F)
  if(length(K) > 1){
    igraph::V(g)$name <- as.character(K)
  }

  igraph::add_edges(g, edges = edges)
}
