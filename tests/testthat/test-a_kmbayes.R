library(testthat)
library(aBKMR)

test_that("sam_py_r and a_kmbayes work correctly", {


  set.seed(111)


  dat <- SimData(n = 500, M = 4)
  y <- dat$y
  Z <- dat$Z
  X <- dat$X


  knots <- sam_py_r(R = Z , nd = as.integer(50), num_nn = as.integer(100), w = FALSE)



  fitkm <- a_kmbayes(y = y, Z = Z, X = X, knots = Z[knots,], iter = 100, verbose = FALSE, varsel = TRUE)


  expect_true(is.list(fitkm))


})
