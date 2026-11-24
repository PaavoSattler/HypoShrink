
test_that("valid companion representations are accepted", {
  H <- rbind(c(1, -1), c(2, -2))
  result <- CompanionHypothesis(H)

  check <- checkCompanion(H, result$L, ytilde = result$ytilde)

  expect_true(check$valid)
  expect_true(check$minimal)
  expect_equal(check$rank, 1)
})

test_that("noncentral companion representations are accepted", {
  H <- rbind(c(1, 0), c(2, 0))
  y <- c(1, 2)
  result <- CompanionHypothesis(H, y)

  check <- checkCompanion(H, result$L, y, result$ytilde)

  expect_true(check$valid)
  expect_true(check$minimal)
})

test_that("incorrect companion representations are rejected", {
  H <- rbind(c(1, 0), c(2, 0))
  y <- c(1, 2)
  result <- CompanionHypothesis(H, y)

  check <- checkCompanion(H, 2 * result$L, y, result$ytilde)

  expect_false(check$valid)
  expect_true(check$minimal)
})

test_that("rank-zero companion representations are handled", {
  H <- matrix(0, 2, 3)
  result <- CompanionHypothesis(H)

  check <- checkCompanion(H, result$L, ytilde = result$ytilde)

  expect_true(check$valid)
  expect_true(check$minimal)
  expect_equal(check$rank, 0)
})

test_that("small matrices are compared relatively", {
  H <- matrix(1e-10, 1, 1)
  L <- 2 * H

  check <- checkCompanion(H, L)

  expect_false(check$valid)
  expect_equal(check$gram_error, 0.75)
})
