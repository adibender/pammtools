context("Predict functions")


 test_that("predict functions work correctly", {

  data("tumor")
  ped <- as_ped(Surv(days, status)~ complications, data = tumor[1:20, ])
  pam <- pamm(ped_status ~ s(tend, k = 3) + complications, data = ped)
  pam2 <- pamm(ped_status ~ s(tend, k = 3) + complications, data = ped,
    engine = "bam", method = "fREML", discrete = TRUE)

  ## predictSurvProb (pec) generic
  spmat <- predictSurvProb.pamm(pam, tumor[21:23,], times = c(90, 500, 1217))
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

  expect_error(predictSurvProb.pamm(pam, tumor[21:23,], times = c(90, 500, 2000)))
  spmat2 <- predictSurvProb.pamm(pam2, tumor[21:23,], times = c(90, 500, 1217))
  expect_identical(round(spmat, 2), round(spmat2, 2))

  }
)
