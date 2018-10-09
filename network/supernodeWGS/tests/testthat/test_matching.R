context("Test matching")

test_that("matching works", {
  set.seed(10)
  name1 <- letters[sample(1:26)]
  name2 <- letters[sample(1:26)]

  res <- matching(name1, name2)

  expect_true(length(res) == 26)
  expect_true(all(res > 0))
  expect_true(all(res <= 26))
  expect_true(all(res %% 1 == 0))
  expect_true(length(res) == length(unique(res)))
})

test_that("matching is correct", {
  set.seed(10)
  name1 <- letters[sample(1:26)]
  name2 <- c(letters[sample(1:26)], as.character(sample(1:1000)))
  name2 <- name2[sample(length(name2))]

  res <- matching(name1, name2)

  expect_true(all(name1 == name2[res]))
})
