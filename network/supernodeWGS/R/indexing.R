#' Pair to index
#'
#' @param idx1 numeric
#' @param idx2 numeric
#' @param max_val numeric
#'
#' @return numeric
#' @export
pair_to_index <- function(idx1, idx2, max_val = 200){
  stopifnot(length(idx1) == 1, length(idx2) == 1)
  stopifnot(idx1 > 0, idx2 > 0, idx1 %% 1 == 0, idx2 %% 1 == 0)
  stopifnot(max_val > 0, max_val %% 1 == 0)
  stopifnot(max_val >= idx1, max_val >= idx2)

  i <- min(c(idx1, idx2)); j <- max(c(idx1, idx2))

  tmp <- ifelse(i > 1, sum(c(max_val:(max_val-(i-2)))), 0)
  tmp + (j-i+1)
}

#' Index to pair
#'
#' @param val numeric
#' @param max_val numeric
#'
#' @return numeric
#' @export
index_to_pair <- function(val, max_val = 200){
  stopifnot(length(val) == 1, val %% 1 == 0, val > 0, val <= max_val*(max_val-1)/2+max_val)

  vec <- c(0, cumsum(c(max_val:1)))
  idx1 <- max(which(vec < val))
  idx2 <- val - vec[idx1] - 1 + idx1
  c(idx1, idx2)
}
