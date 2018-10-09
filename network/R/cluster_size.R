#' Cluster size
#'
#' @param vec vector of test statistics
#' @param clustering vector of clustering assignments
#' @param flag_vec boolean vector of same length as \code{vec}, where indices with
#' \code{TRUE} will be omitted from calculation
#'
#' @return vector
#' @export
cluster_size <- function(vec, clustering, flag_vec){
  stopifnot(length(flag_vec) == length(vec))

  sapply(1:max(clustering), function(x){
    idx <- which(clustering == x)
    length(intersect(which(!is.na(vec[idx])), which(!flag_vec[idx])))
  })
}
