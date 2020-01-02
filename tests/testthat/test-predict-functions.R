context("Predict functions")


 test_that("predict functions work correctly", {

  data("tumor")
  ped <- as_ped(Surv(days, status)~ complications, data = tumor[1:20, ])
  pam <- pamm(ped_status ~ s(tend, k = 3) + complications, data = ped)

  ## predict S3 generic
  predh <- predict(pam, newdata = tumor[21:23, ], type = "hazard")
  predch <- predict(pam, newdata = tumor[21:23, ], type = "cumu_hazard")
  predsp <- predict(pam, newdata = tumor[21:23, ], type = "surv_prob")
  expect_identical(round(predh * 100, 2), c(.19, .02, .02))
  expect_identical(round(predch, 2), c(.2, .96, .96))
  expect_identical(round(predsp, 2), c(.81, .38, .38))

  ## predictSurvProb (pec) generic
  spmat <- predictSurvProb(pam, tumor[21:23,], times = c(90, 500, 1217))
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
  expect_identical(spmat[1,1], predsp[1])
  expect_error(predictSurvProb(pam, tumor[21:23,], times = c(90, 500, 2000)))

  }
)
