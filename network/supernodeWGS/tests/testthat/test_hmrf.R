context("Test hmrf")

## hmrf is correct

test_that("hmrf works", {
  set.seed(10)
  d <- 12
  pval <- stats::runif(d)
  z <- stats::qnorm(1 - pval)
  adj <- as.matrix(huge::huge.generator(n = 50, d = d, graph = "scale-free", verbose = F)$theta)
  seedindex <- rep(0, d)
  seedindex[order(pval, decreasing = F)[1:2]] <- 1

  res <- hmrf(z, adj, seedindex, pthres = 0.5, iter = 2)

  expect_true(is.list(res))
  expect_true(all(c("Iupdate", "post", "b", "c", "mu1", "mu2", "sigma1", "sigma2") %in% names(res)))
})

test_that("hmrf can get c to be large", {
  set.seed(10)
  d <- 12
  adj <- as.matrix(huge::huge.generator(n = 50, d = d, graph = "hub", g = 3, verbose = F)$theta)

  #set one hub to have small p-values, others to have large
  seedindex <- c(1,rep(0, d-1))
  pval <- c(stats::runif(floor(d/3), max = 0.05), stats::runif(d - floor(d/3)))
  z <- stats::qnorm(1 - pval)
  res <- hmrf(z, adj, seedindex, pthres = 0.05, iter = 10)

  #do another case where pval is uninformative
  pval2 <- c(0.01, stats::runif(d-2), 0.01)
  z2 <- stats::qnorm(1 - pval2)
  res2 <- hmrf(z2, adj, seedindex, pthres = 0.05, iter = 10)

  expect_true(res$c >= res2$c)
})

test_that("hmrf gives reasonable posteriors", {
  set.seed(10)
  d <- 12
  adj <- as.matrix(huge::huge.generator(n = 50, d = d, graph = "hub", g = 3, verbose = F)$theta)

  #set one hub to have small p-values, others to have large
  seedindex <- c(1,rep(0, d-1))
  pval <- c(stats::runif(floor(d/3), max = 0.05), stats::runif(d - floor(d/3)))
  z <- stats::qnorm(1 - pval)
  res <- hmrf(z, adj, seedindex, pthres = 0.05, iter = 10)

  expect_true(mean(res$post[1:floor(d/3)]) >= mean(res$post[(floor(d/3)+1):d]))
})
