#' Report results
#'
#' @param vec vector
#' @param posterior vector of numerics
#' @param pvalue vector of numerics
#' @param Iupdate vector of booleans
#'
#' @return data frame with 4 columns
#' @export
report_results <- function(vec, posterior, pvalue, Iupdate){
  stopifnot(length(vec) == length(posterior), length(vec) == length(pvalue),
            length(vec) == length(Iupdate))

  d <- length(posterior)
  rankpost <- sort(posterior)
  localfdr <- sapply(1:d, function(x){mean(rankpost[1:x])})

  flocalfdr <- rep(0,d)
  rankp <- rank(posterior, ties.method="random")
  flocalfdr <- localfdr[rankp] #undo the sorting

  res = data.frame(vec, pvalue, flocalfdr, Iupdate)
  names(res) = c("Name","p.value","FDR", "indicator")

  res
}
