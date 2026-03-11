context("Predict functions")


test_that("predict functions work correctly", {
  data("tumor")
  ped <- as_ped(Surv(days, status) ~ complications, data = tumor[1:20, ])
  pam <- pamm(ped_status ~ s(tend, k = 3) + complications, data = ped)
  pam2 <- pamm(
    ped_status ~ s(tend, k = 3) + complications,
    data = ped,
    engine = "bam",
    method = "fREML",
    discrete = TRUE
  )

  ## predictSurvProb (pec) generic
  spmat <- predictSurvProb.pamm(pam, tumor[21:23, ], times = c(90, 500, 1217))
  expect_identical(
    round(spmat, 2),
    matrix(
      c(
        rep(.81, 3),
        rep(.46, 3),
        rep(.38, 3)
      ),
      nrow = 3,
      ncol = 3
    )
  )

  expect_error(predictSurvProb.pamm(
    pam,
    tumor[21:23, ],
    times = c(90, 500, 2000)
  ))
  spmat2 <- predictSurvProb.pamm(pam2, tumor[21:23, ], times = c(90, 500, 1217))
  expect_identical(round(spmat, 2), round(spmat2, 2))
})

test_that("get_terms returns tidy partial effects for requested terms", {
  z95 <- stats::qnorm(0.975)
  fit <- survival::coxph(
    survival::Surv(time, status) ~
      survival::pspline(karno) + survival::pspline(age),
    data = survival::veteran
  )

  terms_df <- get_terms(
    data = survival::veteran,
    fit = fit,
    terms = c("karno", "age")
  )

  expect_true(is.data.frame(terms_df))
  expect_identical(
    names(terms_df),
    c("term", "x", "eff", "se", "ci_lower", "ci_upper")
  )
  expect_equal(sort(unique(terms_df$term)), c("age", "karno"))
  expect_equal(as.integer(table(terms_df$term)[["age"]]), 100L)
  expect_equal(as.integer(table(terms_df$term)[["karno"]]), 100L)
  expect_true(all(terms_df$ci_lower <= terms_df$eff))
  expect_true(all(terms_df$ci_upper >= terms_df$eff))
  expect_equal(
    terms_df$ci_upper - terms_df$eff,
    z95 * terms_df$se
  )
})
