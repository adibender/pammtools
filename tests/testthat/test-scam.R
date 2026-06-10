context("Support for shape-constrained additive models (scam)")

library(scam)

data("tumor")
ped <- tumor[1:200, ] %>%
  as_ped(
    Surv(days, status) ~ age + complications,
    cut = c(0, 50, 100, 200, 300, 400)
  )

mpam <- scam(
  ped_status ~ s(tend, k = 5, bs = "mpd") + complications,
  data = ped,
  family = poisson(),
  offset = ped$offset
)
mpam2 <- scam(
  ped_status ~ s(tend, k = 5, bs = "mpd") + complications + age,
  data = ped,
  family = poisson(),
  offset = ped$offset
)

test_that("hazards and standard errors match predict.scam", {
  nd <- ped_info(ped)
  haz <- add_hazard(nd, mpam)
  eta <- as.numeric(predict(mpam, newdata = nd, type = "link"))
  expect_equal(haz$hazard, exp(eta), tolerance = 1e-8)
  se <- as.numeric(
    predict(mpam, newdata = nd, type = "link", se.fit = TRUE)$se.fit
  )
  expect_equal(haz$se, se, tolerance = 1e-6)
})

test_that("monotonicity constraint is respected", {
  haz <- add_hazard(ped_info(ped), mpam)
  expect_true(all(diff(haz$hazard) <= 1e-10))
})

test_that("hazard functions work for scam", {
  expect_data_frame(
    haz <- add_hazard(ped_info(ped), mpam),
    nrows = 5L,
    ncols = 11L
  )
  expect_true(all(haz$ci_lower <= haz$hazard))
  expect_true(all(haz$ci_upper >= haz$hazard))
  expect_error(add_hazard(haz, mpam))
  expect_data_frame(
    add_hazard(haz, mpam, overwrite = TRUE),
    nrows = 5L,
    ncols = 11L
  )

  haz2 <- add_hazard(ped_info(ped), mpam, type = "link")
  expect_true(all(haz2$ci_lower < haz2$hazard))
  expect_true(all(haz2$ci_upper > haz2$hazard))

  ## delta rule
  haz3 <- add_hazard(ped_info(ped), mpam, ci_type = "delta")
  expect_data_frame(haz3, nrows = 5L, ncols = 11L)
  expect_true(all(haz3$ci_lower <= haz3$hazard))
  expect_true(all(haz3$ci_upper >= haz3$hazard))

  ## simulation based ci
  set.seed(123)
  haz4 <- add_hazard(ped_info(ped), mpam, ci_type = "sim")
  expect_data_frame(haz4, nrows = 5L, ncols = 11L)
  expect_true(all(haz4$ci_lower <= haz4$hazard))
  expect_true(all(haz4$ci_upper >= haz4$hazard))

  ## hazard with reference (i.e. hazard ratio)
  hr <- add_hazard(ped_info(ped), mpam2, reference = list(age = c(30)))
  expect_data_frame(hr, nrows = 5L, ncols = 11L)
  # hr = 1 if reference = data
  hr2 <- ped_info(ped) %>%
    add_hazard(mpam2, reference = list(age = mean(.$age)))
  expect_equal(hr2$hazard, rep(1, 5))
})

test_that("cumulative hazard and survival probabilities work for scam", {
  for (ci_type in c("default", "delta", "sim")) {
    set.seed(123)
    cumu <- add_cumu_hazard(ped_info(ped), mpam, ci_type = ci_type)
    expect_subset(
      c("cumu_hazard", "cumu_lower", "cumu_upper"),
      colnames(cumu)
    )
    expect_true(all(diff(cumu$cumu_hazard) >= 0))
    expect_true(all(cumu$cumu_lower <= cumu$cumu_hazard))
    expect_true(all(cumu$cumu_upper >= cumu$cumu_hazard))

    set.seed(123)
    surv <- add_surv_prob(ped_info(ped), mpam, ci_type = ci_type)
    expect_subset(c("surv_prob", "surv_lower", "surv_upper"), colnames(surv))
    expect_true(all(surv$surv_prob >= 0 & surv$surv_prob <= 1))
    expect_true(all(diff(surv$surv_prob) <= 0))
    expect_true(all(surv$surv_lower <= surv$surv_prob))
    expect_true(all(surv$surv_upper >= surv$surv_prob))
  }
})

test_that("grouped cumulative hazards work for scam", {
  cumu <- ped %>%
    make_newdata(tend = unique(tend), complications = unique(complications)) %>%
    group_by(complications) %>%
    add_cumu_hazard(mpam)
  expect_data_frame(cumu, nrows = 10L)
  cumu_split <- split(cumu$cumu_hazard, cumu$complications)
  expect_true(all(vapply(cumu_split, function(x) all(diff(x) >= 0), logical(1))))
})

test_that("add_term works for scam", {
  term <- ped %>%
    make_newdata(age = seq_range(age, 10)) %>%
    add_term(mpam2, term = "age", reference = list(age = mean(.$age)))
  expect_data_frame(term, nrows = 10L)
  expect_subset(c("fit", "ci_lower", "ci_upper"), colnames(term))
  expect_true(all(term$ci_lower <= term$fit))
  expect_true(all(term$ci_upper >= term$fit))
})

test_that("pamm works with engine scam", {
  pam_scam <- pamm(
    ped_status ~ s(tend, k = 5, bs = "mpd") + complications,
    data = ped,
    engine = "scam"
  )
  expect_class(pam_scam, c("pamm", "scam"))
  haz <- add_hazard(ped_info(ped), pam_scam)
  expect_true(all(diff(haz$hazard) <= 1e-10))
  expect_equal(haz$hazard, add_hazard(ped_info(ped), mpam)$hazard,
    tolerance = 1e-6)

  # compatibility with pec::predictSurvProb
  surv_prob <- pec::predictSurvProb(
    pam_scam,
    newdata = tumor[3:6, ],
    times = c(100, 300)
  )
  expect_matrix(surv_prob, nrows = 4L, ncols = 2L)
  expect_true(all(surv_prob >= 0 & surv_prob <= 1))
  expect_true(all(surv_prob[, 2] <= surv_prob[, 1]))
})

test_that("tidiers work for scam", {
  fixed <- tidy_fixed(mpam2)
  expect_data_frame(fixed, nrows = 2L, ncols = 4L)
  expect_identical(
    colnames(fixed),
    c("variable", "coef", "ci_lower", "ci_upper")
  )

  smooth <- tidy_smooth(mpam)
  expect_data_frame(smooth, nrows = 100L)
  expect_subset(c("x", "fit", "ci_lower", "ci_upper"), colnames(smooth))
})

test_that("cumulative coefficients work for scam", {
  set.seed(123)
  cumu_coef <- get_cumu_coef(mpam2, ped, terms = c("age", "complications"))
  expect_data_frame(cumu_coef, nrows = 10L, ncols = 6L)
  expect_subset(
    c("method", "variable", "time", "cumu_hazard", "cumu_lower", "cumu_upper"),
    colnames(cumu_coef)
  )
})

test_that("gg_smooth works for scam", {
  gg <- gg_smooth(ped, mpam, terms = "tend")
  expect_is(gg, "ggplot")
})

test_that("CIF works for scam", {
  set.seed(211758)
  df <- data.frame(
    time = rexp(20),
    status = sample(c(0, 1, 2), 20, replace = TRUE)
  )
  ped_cr <- as_ped(df, Surv(time, status) ~ ., id = "id") %>%
    mutate(cause = as.factor(cause))
  scam_cr <- scam(
    ped_status ~ s(tend, k = 5, bs = "mpd", by = cause),
    data = ped_cr,
    family = poisson(),
    offset = ped_cr$offset
  )
  set.seed(123)
  ndf <- ped_cr %>%
    make_newdata(tend = unique(tend), cause = unique(cause)) %>%
    group_by(cause) %>%
    add_cif(scam_cr, nsim = 50L)
  expect_subset(c("cif", "cif_lower", "cif_upper"), colnames(ndf))
  expect_true(all(ndf$cif >= ndf$cif_lower & ndf$cif <= ndf$cif_upper))
  expect_true(all(ndf$cif >= 0 & ndf$cif <= 1))
})
