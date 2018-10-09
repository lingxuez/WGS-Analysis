context("Test for correlation")

## cor_test is correct

test_that("cor_test works", {
  set.seed(10)
  res <- cor_test(runif(10))

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res >= 0)
  expect_true(res <= 1)
})

test_that("cor_test gives relatively correct p-values", {
  trials <- 100
  vec1 <- sapply(1:trials, function(x){
    set.seed(x)
    mat1 <- MASS::mvrnorm(100, rep(0, 20), diag(20))
    mat2 <- MASS::mvrnorm(100, rep(0, 20), diag(20))
    cor_mat <- stats::cor(mat1, mat2)
    cor_test(as.numeric(cor_mat))
  })

  vec2 <- sapply(1:trials, function(x){
    set.seed(x)
    mat1 <- MASS::mvrnorm(100, rep(0, 20), diag(20))
    mat2 <- mat1 + 0.1*MASS::mvrnorm(100, rep(0, 20), diag(20))
    cor_mat <- stats::cor(mat1, mat2)
    cor_test(as.numeric(cor_mat))
  })

  expect_true(sum(abs(stats::quantile(vec1) - seq(0,1,length.out=5))) <=
                sum(abs(stats::quantile(vec2) - seq(0,1,length.out=5))))
})
