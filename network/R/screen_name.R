#' Screen names based on words in certain positions
#'
#' @param vec vector of names
#' @param position index
#' @param stopwords words to remove
#'
#' @return vector of booleans
#' @export
screen_name <- function(vec, position = 4,
                        stopwords = c("Any", "CodingRegion", "FrameshiftRegion", "LoFRegion",
                                  "MissenseRegion", "SilentRegion", "MissenseHVARDRegionSimple")){
  vec <- as.character(vec)

  bool_vec <- sapply(vec, function(x){
    words <- strsplit(x, split = "_")[[1]]
    !any(stopwords %in% words[position])
  })

  as.vector(bool_vec)
}

#' Tabulate names based on position
#'
#' @param vec vector of names
#' @param position index
#'
#' @return vector of counts
#' @export
tabulate_names <- function(vec, position = 4){
  vec <- as.character(vec)

  tab <- sapply(vec, function(x){
    words <- strsplit(x, split = "_")[[1]]
    words[position]
  })

  table(tab)
}
