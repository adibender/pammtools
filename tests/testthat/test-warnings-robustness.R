context("Warnings robustness")

test_that("direct warning methods enforce interval mismatch guardrails", {
  ped <- get_tumor_ped_fixture()

  pem <- glm(
    ped_status ~ 0 + interval + complications,
    data = ped,
    family = poisson(),
    offset = offset
  )
  expect_error(
    pammtools:::warn_about_new_time_points.glm(
      object = pem,
      newdata = data.frame(interval = factor("(999,1000]")),
      time_var = "interval"
    ),
    "Time points/intervals in new data not equivalent"
  )

  pam <- pamm(ped_status ~ s(tend, k = 5) + complications, data = ped)
  expect_warning(
    pammtools:::warn_about_new_time_points.pamm(
      object = pam,
      newdata = data.frame(interval = factor("(999,1000]"))
    ),
    "Time points/intervals in new data not equivalent"
  )
})

test_that("generic warn_about_new_time_points dispatches for glm/pamm objects", {
  ped <- get_tumor_ped_fixture()

  pem <- glm(
    ped_status ~ 0 + interval + complications,
    data = ped,
    family = poisson(),
    offset = offset
  )
  pam <- pamm(ped_status ~ s(tend, k = 5) + complications, data = ped)

  expect_error(
    pammtools:::warn_about_new_time_points(
      object = pem,
      newdata = data.frame(interval = factor("(999,1000]")),
      time_var = "interval"
    ),
    "Time points/intervals in new data not equivalent"
  )

  expect_warning(
    pammtools:::warn_about_new_time_points(
      object = pam,
      newdata = data.frame(interval = factor("(999,1000]"))
    ),
    "Time points/intervals in new data not equivalent"
  )
})

test_that("as_ped fails fast if no event is present", {
  all_censored <- data.frame(
    id = 1:5,
    time = 1:5,
    status = 0,
    x = c(0.2, 1.1, -0.3, 0.9, 0.4)
  )

  expect_error(
    as_ped(
      all_censored,
      survival::Surv(time, status) ~ x,
      id = "id",
      cut = 0:5
    ),
    "No events in data! Check your status variable."
  )
})
