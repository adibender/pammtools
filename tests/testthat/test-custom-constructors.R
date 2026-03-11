context("Custom smooth constructors")

test_that("fdl constructor sets penalty blocks and ranks consistently", {
  base <- construct_fdl_smooth(xt = list())
  fullrank <- construct_fdl_smooth(xt = list(fullrankpen = TRUE))
  lead <- construct_fdl_smooth(xt = list(ridge = TRUE, leadpen = TRUE))
  lag <- construct_fdl_smooth(xt = list(ridge = TRUE, lagpen = TRUE))
  both <- construct_fdl_smooth(
    xt = list(ridge = TRUE, leadpen = TRUE, lagpen = TRUE)
  )
  constrained <- construct_fdl_smooth(xt = list(constrain = TRUE))

  expect_equal(base$rank, 6)
  expect_equal(fullrank$rank, 8)
  expect_equal(lead$rank, c(6, 4))
  expect_equal(lag$rank, c(6, 4))
  expect_equal(both$rank, c(6, 8))

  expect_equal(which(diag(lead$S[[2]]) > 0), 1:4)
  expect_equal(which(diag(lag$S[[2]]) > 0), 5:8)
  expect_equal(which(diag(both$S[[2]]) > 0), 1:8)

  expect_true(!is.null(constrained$C))
  expect_true(isTRUE(attr(constrained$C, "always.apply")))
  expect_equal(
    as.numeric(constrained$C[1, ]),
    as.numeric(constrained$X[1, ]),
    tolerance = 1e-12
  )
})

test_that("fdl basis is discoverable by mgcv smooth constructor", {
  set.seed(611)
  dat <- data.frame(
    y = rpois(100, lambda = 3),
    x = runif(100)
  )

  expect_no_error(
    mgcv::gam(
      y ~ s(x, bs = "fdl", k = 8),
      family = poisson(),
      data = dat
    )
  )
})
