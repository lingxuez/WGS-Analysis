#' Matching two vectors
#'
#' Returns indicies so \code{all(name1 == name2[matching(name1, name2)])}.
#'
#' @param name1 vector
#' @param name2 vector
#' @param safety boolean
#'
#' @return numeric vector
#' @export
matching <- function(name1, name2, safety = T){
  if(safety) stopifnot(all(name1 %in% name2))

  idx <- which(name1 %in% name2)

  vec <- rep(NA, length(name1))
  vec[idx] <- as.numeric(plyr::mapvalues(name1[idx], from = name2, to = 1:length(name2),
                         warn_missing = F))

  vec
}

#' Approximate matching function
#'
#' This function creates a vector that is length of \code{name1} where the
#' \code{i}th element of this vector contains the value (say, called \code{x}) such
#' that \code{name1[i]} is matched to \code{name2[x]}.
#'
#' There are few tricks this function does. First, it splits all strings in
#' \code{name1} and \code{name2} by underscores, (i.e., "_"). Next, it removes
#' all words with less than 3 characters. Then, all subwords are removed.
#' (For example, if the two words in the string are "Antisense" and "AntisenseMutation",
#' the word "Antisense" is removed.)
#'
#' Assuming there is only one string in both \code{name1} and \code{name2} after
#' this processing with no words (i.e., empty string), the empty strings are matched
#' to each other.
#'
#' Finally, the element in \code{name2} is returned
#' such 1) this element has the highest percentage
#' of its words in \code{name1[i]} of all elements in \code{name2}
#' and 2) \code{name1[i]} has highest percentage of its words in this element
#' among all elements in \code{name2}.
#' If no elements satisfy this condition (i.e., the intersection), then
#' an \code{NA} is returned. Otherwise, if there are more than one element
#' that satifies this condition, a \code{0} is returned.
#'
#' @param name1 vector
#' @param name2 vector
#' @param verbose boolean
#'
#' @return numeric vector
#' @export
approximate_matching <- function(name1, name2, verbose = T){
  transformation <- function(x){
    vec <- sort(strsplit(x, split = "_")[[1]])
    vec <- vec[which(sapply(vec, nchar) > 3)]
    if(length(vec) > 0) vec <- .remove_subwords(vec)
    sort(vec)
  }

  name1_transform <- lapply(name1, transformation)
  name2_transform <- lapply(name2, transformation)

  name1_vec <- as.vector(sapply(name1_transform, paste0, collapse = "_"))
  name2_vec <- as.vector(sapply(name2_transform, paste0, collapse = "_"))

  stopifnot(length(name1_vec) == length(name1_transform))
  stopifnot(length(name2_vec) == length(name2_transform))

  match_vec <- sapply(1:length(name1_transform), function(x){
    if(verbose & x %% floor(length(name1_transform)/10) == 0) cat('*')

    if(name1_vec[x] %in% name2_vec) {
      idx <- which(name2_vec == name1_vec[x])
      stopifnot(length(idx) == 1)
      return(idx)
    }

    if(length(name1_transform[[x]]) == 0) {
      idx <- which(sapply(name2_transform, length) == 0)
      stopifnot(length(idx) == 1)
      return(idx)
    }

    vec <- sapply(name2_transform, function(y){
      leny <- ifelse(length(y) != 0, length(y), 1)
      c(length(which(name1_transform[[x]] %in% y))/length(name1_transform[[x]]),
        length(which(y %in% name1_transform[[x]]))/leny)
    })

    idx1 <- which(vec[1,] == max(vec[1,]))
    idx2 <- which(vec[2,] == max(vec[2,]))
    idx <- intersect(idx1, idx2)

    if(length(idx) == 0) return(NA)
    ifelse(length(idx) == 1, idx, 0)
  })

  match_vec
}

.remove_subwords <- function(vec){
  len <- length(vec)
  bool_vec <- sapply(1:len, function(x){
    any(grepl(vec[x], vec[-x]))
  })

  vec[which(!bool_vec)]
}
