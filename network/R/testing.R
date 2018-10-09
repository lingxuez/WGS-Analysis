#' Correlation test
#'
#' @param vec numeric
#'
#' @return scalar
#' @export
cor_test <- function(vec){
  n <- length(vec)
  vec <- atanh(vec)

  se <- stats::sd(vec)/sqrt(n)
  z <- abs(mean(vec)/se)
  2*stats::pnorm(-z)
}
