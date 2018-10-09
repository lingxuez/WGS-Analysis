context("Test indexing")

## pair_to_index is correct

test_that("pair_to_index works", {
  res <- pair_to_index(4,2)
  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 1)
  expect_true(res > 0)
  expect_true(res <= 200*199/2+200)
  expect_true(res %% 1 == 0)
})

test_that("pair_to_index gives unique values", {
  max_val <- 20
  combn_mat <- combn(max_val,2)
  combn_mat <- cbind(combn_mat, sapply(1:max_val, function(x){rep(x,2)}))
  res <- apply(combn_mat, 2, function(x){
    pair_to_index(x[1],x[2],max_val)
  })

  expect_true(length(res) == length(unique(res)))
  expect_true(all(sort(res) == 1:(ncol(combn_mat))))
})

#####################

## index_to_pair is correct

test_that("index_to_pair works", {
  res <- index_to_pair(132, max_val = 20)

  expect_true(is.numeric(res))
  expect_true(!is.matrix(res))
  expect_true(length(res) == 2)
  expect_true(all(res > 0))
  expect_true(all(res %% 1 == 0))
})

test_that("index_to_pair inverses pair_to_index", {
  max_val <- 20
  combn_mat <- combn(max_val,2)
  combn_mat <- cbind(combn_mat, sapply(1:max_val, function(x){rep(x,2)}))
  vec <- apply(combn_mat, 2, function(x){
    pair_to_index(x[1],x[2],max_val)
  })

  res <- sapply(vec, function(x){index_to_pair(x, max_val)})

  expect_true(all(res == combn_mat))
})
