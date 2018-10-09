#partial likelihood function
.partial_likelihood <- function(b, c, graph_term, i_vec){
  d <- length(i_vec)

  fv.fun <- function(i){
    new1 <- exp(b*i_vec[i]+c*i_vec[i]*graph_term[i])
    new2 <- exp(b*(1-i_vec[i])+c*(1-i_vec[i])*graph_term[i])
    fvalue <- log(new1/(new1+new2))
  }

  sum(sapply(1:d,fv.fun))
}

##estimate Ising model parameters b,c
.optimize_bc <- function(func, adj, i_vec, times, tol = 1e-5, range = c(-10, 10)){
  b <- 0; c <- 1

  graph_term <- i_vec%*%adj
  for (k in 1:times){
    b_new <- stats::optimize(func, range, c = c, graph_term = graph_term,
                             i_vec = i_vec, maximum = T)$maximum
    c_new <- stats::optimize(func, range, b = b_new, graph_term = graph_term,
                             i_vec = i_vec, maximum = T)$maximum
    if (abs(b_new-b) < tol & abs(c_new-c) < tol){
      break()
    }

    b <- b_new; c <- c_new
  }

  list(b = b, c = c)
}


#' Hidden markov random field
#'
#' Use a Ising model to update the pvalues according to the graphical structure.
#' There are two groups, indices in group 0 have mean 0. Indices in group 1 have
#' mean greater than 0.
#'
#' In the output, \code{post} refers to the posterior probability of being in
#' group 1.
#'
#' @param z a vector of d z-values
#' @param adj a dxd adjancey matrix (0,1 entries)
#' @param seedindex a (0,1) vector of length d, where 1 means the gene is in group 1.
#' @param null_mean the mean of the null distribution
#' @param null_sigma the sd of the null distribution
#' @param pthres threshold for p-values to serve as an initialization
#' @param iter number of iterations
#' @param verbose boolean
#' @param tol tolerance for values that should be treated as zero in absolute value
#'
#' @return a list
#' @export
hmrf <- function(z, adj, seedindex = rep(0, length(z)),
                 null_mean = 0, null_sigma = NA,
                 pthres = 0.05, iter = 100,
                 verbose = FALSE, tol = 1e-3){
  stopifnot(length(z) == length(seedindex), all(dim(adj) == length(z)))

  #permute all the entries to avoid emphasis on current order
  d <- length(z)
  idx <- sample(1:d)
  z <- z[idx]; adj <- adj[idx,idx]; seedindex <- seedindex[idx]

  i_vec <- as.numeric(z > stats::qnorm(1 - pthres, mean = null_mean, sd = ifelse(is.na(null_sigma), 1, null_sigma)))
  b <- 0; c <- 0

  mu2 <- mean(z[i_vec==1 & seedindex==0])
  mu1 <- null_mean
  if(is.na(null_sigma)){
    sigma2 <- stats::var(z)
    sigma1 <- sigma2
  } else {
    sigma2 <- null_sigma
    sigma1 <- null_sigma
  }
  posterior <- rep(0,d)

  for (iteri in 1:iter){

    if(verbose) print(paste("Start iteration: ", iteri))
    res <- .optimize_bc(.partial_likelihood, adj, i_vec, 20)
    b_new <- res$b; c_new <- res$c

    if (abs(c-c_new)<tol & abs(b-b_new)<tol)  break()
    b <- b_new; c <- c_new

    for (i in 1:d){
      new1 <- exp(b*i_vec[i] + c*i_vec[i]*adj[i,]%*%i_vec)
      new2 <- exp(b*(1-i_vec[i]) + c*(1-i_vec[i])*adj[i,]%*%i_vec)
      p1 <- stats::dnorm(z[i], mu2*i_vec[i] + mu1*(1-i_vec[i]), sqrt(sigma2*i_vec[i] + sigma1*(1-i_vec[i]))) * new1/(new1+new2)
      p2 <- stats::dnorm(z[i], mu2*(1-i_vec[i]) + mu1*(i_vec[i]), sqrt(sigma2*(1-i_vec[i]) + sigma1*(i_vec[i]))) * new2/(new1+new2)

      if (i_vec[i] == 1){
        posterior[i] <- p1/(p1+p2)
      } else {
        posterior[i] <- p2/(p1+p2)
      }

      if (p2 > p1){ i_vec[i] <- 1-i_vec[i] }
      if (seedindex[i] != 0){ i_vec[i] <- 1 }
    }

    mu2 <- sum(posterior[seedindex==0]*z[seedindex==0])/sum(posterior[seedindex==0])
    sigma2 <- sum(posterior[seedindex==0]*(z[seedindex==0]-mu2)^2)/sum(posterior[seedindex==0])
    if(is.na(null_sigma)) {
      sigma1 <- sum((1-posterior[seedindex==0])*(z[seedindex==0])^2)/sum(1-posterior[seedindex==0])
      sigmas <- (sigma1*sum(posterior[seedindex==0])+sigma2*sum(1-posterior[seedindex==0])) / length(posterior)
      sigma2 <- sigmas; sigma1 <- sigmas
    }

    if(verbose) {
      print(paste0("Iteration: ",iteri," has ", sum(i_vec)," genes set with Iupdate = 1."))
      print(paste0("Parameters: ", round(mu1,3), ", ", round(mu2,3), "//", round(sigma1,3), ",", round(sigma2,3)))
    }
  }

  #un-permute the elements
  i_vec <- i_vec[order(idx)]; posterior <- posterior[order(idx)]

  list(Iupdate = i_vec, post = posterior, b = b, c = c, mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2)
}
