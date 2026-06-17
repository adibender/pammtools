context("Interval-censored data via multiple imputation")

make_ic_data <- function(n = 200, seed = 1, cut = seq(0, 10, by = 0.25)) {
  set.seed(seed)
  df <- data.frame(x = runif(n, -1, 1))
  sdf <- sim_pexp(~ -2 + 0.5 * x, df, cut = cut)
  add_inspections(sdf, rate = 1, max_time = max(cut))
}

test_that("interval-censored Surv is detected and parsed", {
  icd <- make_ic_data(120)
  f <- Surv(L, R, type = "interval2") ~ x
  expect_identical(detect_ic(f, icd), "interval2")
  # right-censored / left-truncated responses are NOT treated as IC
  expect_identical(
    detect_ic(Surv(true_time, rep(1, nrow(icd))) ~ x, icd),
    "none"
  )

  ic <- parse_ic_surv(f, icd)
  expect_true(all(c("ic_L", "ic_R", "ic_kind", "id") %in% names(ic)))
  expect_true(all(
    as.character(ic$ic_kind) %in%
      c("exact", "right", "left", "interval")
  ))
  # bounds bracket the truth
  ev <- as.character(ic$ic_kind) %in% c("interval", "left")
  expect_true(all(icd$true_time[ev] > ic$ic_L[ev] - 1e-9))
  expect_true(all(
    icd$true_time[ev] <= ic$ic_R[ev] + 1e-9 |
      is.infinite(ic$ic_R[ev])
  ))
})

test_that("as_ped dispatches interval-censored data to the IC pipeline", {
  icd <- make_ic_data(120)
  cut <- seq(0, 10, by = 0.5)
  ped <- as_ped(icd, Surv(L, R, type = "interval2") ~ x, cut = cut)
  expect_s3_class(ped, "ped_ic_init")
  expect_s3_class(ped, "ped")
  expect_identical(attr(ped, "breaks"), cut)
  expect_true(!is.null(attr(ped, "ic")))
  # standard right-censored data is unaffected
  std <- as_ped(tumor[1:50, ], Surv(days, status) ~ age)
  expect_false(inherits(std, "ped_ic_init"))
})

test_that("IC cuts must cover finite event upper bounds", {
  df <- data.frame(
    x = c(0, 1, 0),
    L = c(1, 2, 5),
    R = c(2, 6, Inf)
  )
  expect_error(
    as_ped_ic(df, Surv(L, R, type = "interval2") ~ x, cut = 0:5),
    "cover all finite event upper bounds"
  )
})

test_that("conditional-hazard sampler draws within (L, R]", {
  icd <- make_ic_data(200, seed = 3)
  cut <- seq(0, 10, by = 0.5)
  fit <- pamm_ic(Surv(L, R, type = "interval2") ~ x, icd, cut = cut, m = 1)
  ic <- fit$ic
  cache <- ic_pred_cache(fit$init_fit, ic, cut)
  for (rep in 1:5) {
    ti <- impute_ic_times(fit$init_fit, ic, cut, cache = cache)
    idx <- which(as.character(ic$ic_kind) %in% c("interval", "left"))
    expect_true(all(ti[idx] > ic$ic_L[idx] - 1e-8))
    expect_true(all(ti[idx] <= pmin(ic$ic_R[idx], max(cut)) + 1e-8))
    # exact / right rows are passed through unchanged
    keep <- which(as.character(ic$ic_kind) %in% c("exact", "right"))
    expect_equal(ti[keep], ic$ic_L[keep])
  }
})

test_that("sampler is calibrated (probability integral transform is uniform)", {
  # If T_i is drawn from the conditional CDF F_i on (L_i, R_i], then
  # U_i = F_i(T_i) ~ Uniform(0, 1). Test this across all (heterogeneous)
  # interval-censored subjects with a single KS test.
  icd <- make_ic_data(400, seed = 11)
  cut <- seq(0, 10, by = 0.5)
  fit <- pamm_ic(Surv(L, R, type = "interval2") ~ x, icd, cut = cut, m = 1)
  ic <- fit$ic
  cache <- ic_pred_cache(fit$init_fit, ic, cut)

  ii <- cache$ii
  n_int <- nrow(ii)
  hm <- matrix(as.numeric(exp(cache$X %*% coef(fit$init_fit))), nrow = n_int)
  Hcut <- rbind(0, apply(hm * ii$intlen, 2, cumsum))
  cutv <- c(ii$tstart[1], ii$tend)
  evalH <- function(t, s) {
    j <- pmin(
      pmax(
        findInterval(t, cutv, left.open = TRUE, rightmost.closed = TRUE),
        1L
      ),
      n_int
    )
    Hcut[cbind(j, s)] + hm[cbind(j, s)] * (t - cutv[j])
  }

  ti <- impute_ic_times(fit$init_fit, ic, cut, cache = cache)
  idx <- which(as.character(ic$ic_kind) %in% c("interval", "left"))
  HL <- evalH(ic$ic_L[idx], idx)
  HR <- evalH(pmin(ic$ic_R[idx], max(cut)), idx)
  Ht <- evalH(ti[idx], idx)
  U <- (exp(-HL) - exp(-Ht)) / (exp(-HL) - exp(-HR))
  U <- U[is.finite(U)]
  ks <- suppressWarnings(ks.test(U, "punif"))
  expect_gt(ks$p.value, 0.01)
})

test_that("pamm_ic returns m ordinary pamm fits and prints", {
  icd <- make_ic_data(150, seed = 4)
  fit <- pamm_ic(
    Surv(L, R, type = "interval2") ~ x,
    icd,
    cut = seq(0, 10, by = 0.5),
    m = 5
  )
  expect_s3_class(fit, "pamm_ic")
  expect_length(fit$fits, 5)
  expect_true(all(vapply(fit$fits, inherits, logical(1), "pamm")))
  expect_output(print(fit), "interval-censored")
})

test_that("pamm_ic resolves dot formulas before fitting", {
  icd <- make_ic_data(150, seed = 12)
  icd <- icd[, c("id", "x", "L", "R")]
  fit <- pamm_ic(
    Surv(L, R, type = "interval2") ~ .,
    icd,
    cut = seq(0, 10, by = 0.5),
    m = 1
  )

  expect_s3_class(fit, "pamm_ic")
  expect_false("." %in% all.vars(fit$model_formula))
  expect_true("x" %in% all.vars(fit$model_formula))
})

test_that("pamm_ic stores slimmed fits and a pooled summary", {
  icd <- make_ic_data(250, seed = 4)
  fit <- pamm_ic(
    Surv(L, R, type = "interval2") ~ x,
    icd,
    cut = seq(0, 10, by = 0.5),
    m = 10
  )

  # fits are slimmed: no per-observation slots, but predict/coef/vcov still work
  expect_null(fit$fits[[1]]$model)
  expect_null(fit$fits[[1]]$residuals)
  expect_true(!is.null(coef(fit$fits[[1]])))
  expect_equal(dim(vcov(fit$fits[[1]])), rep(length(coef(fit$fits[[1]])), 2))
  nd <- make_newdata(
    as_ped(icd, Surv(L, R, type = "interval2") ~ x, cut = seq(0, 10, by = 0.5)),
    tend = unique(tend)
  )
  expect_silent(predict(fit$fits[[1]], nd, type = "lpmatrix"))

  # pooled object: summary container, not a fabricated gam fit
  expect_true(!is.null(fit$pooled))
  expect_s3_class(fit$pooled, "pamm_ic_pool")
  expect_false(inherits(fit$pooled, "gam"))
  expect_setequal(
    names(coef(fit$fits[[1]]))[seq_len(fit$fits[[1]]$nsdf)],
    rownames(fit$pooled$p.table)
  )
  expect_true(!is.null(fit$pooled$p.table))
  expect_true(all(c("riv", "fmi", "fmi.table") %in% names(fit$pooled)))
  expect_setequal(rownames(fit$pooled$fmi.table), rownames(fit$pooled$p.table))
  expect_setequal(rownames(fit$pooled$smooth.fmi), rownames(fit$pooled$s.table))
  expect_true(all(is.finite(fit$pooled$smooth.fmi)))
  # pooled (within + between) variance is >= mean within-imputation variance
  within <- mean(vapply(fit$fits, function(f) vcov(f)["x", "x"], 0))
  expect_gte(fit$pooled$Vp["x", "x"], within - 1e-10)
})

test_that("QOI draws receive Rubin finite-imputation inflation", {
  estimates <- list(1, 2, 4)
  draws <- lapply(estimates, function(x) matrix(x, nrow = 1L))

  inflated <- rubin_inflate_qoi_draws(draws, estimates)
  inflated_values <- as.numeric(do.call(cbind, inflated))
  b <- stats::var(unlist(estimates))
  rubin_between <- (1 + 1 / length(estimates)) * b

  expect_equal(
    mean((inflated_values - mean(inflated_values))^2),
    rubin_between,
    tolerance = 1e-12
  )
  expect_equal(rubin_inflate_qoi_draws(draws[1], estimates[1]), draws[1])
})

test_that("print/summary report the pooled summary", {
  icd <- make_ic_data(200, seed = 6)
  fit <- pamm_ic(
    Surv(L, R, type = "interval2") ~ x,
    icd,
    cut = seq(0, 10, by = 0.5),
    m = 5
  )
  expect_output(print(fit), "interval-censored")
  expect_output(print(fit), "via multiple imputation")
  # print() is a compact header only; the coefficient tables live in summary()
  expect_output(print(fit), "Use summary\\(\\) for pooled inference")
  expect_failure(expect_output(print(fit), "Pooled parametric coefficients"))
  s <- summary(fit)
  expect_s3_class(s, "summary.pamm_ic")
  expect_true(!is.null(s$p.table) && !is.null(s$s.table))
  expect_true(!is.null(s$fmi.table) && !is.null(s$smooth.fmi))
  expect_output(print(s), "Rubin")
  expect_output(print(s), "Fraction of missing information")
  expect_output(print(s), "Smooth terms")
})

test_that("exact-only data reproduces a plain right-censored PAMM", {
  set.seed(9)
  df <- data.frame(x = runif(150, -1, 1))
  ex <- sim_pexp(~ -2 + 0.5 * x, df, cut = seq(0, 10, by = 0.5))
  ex$L <- ex$time
  ex$R <- ifelse(ex$status == 1, ex$time, Inf)
  cut <- seq(0, 10, by = 0.5)
  fic <- pamm_ic(Surv(L, R, type = "interval2") ~ x, ex, cut = cut, m = 2)
  ref <- pamm(
    ped_status ~ s(tend) + x,
    data = as_ped(ex, Surv(time, status) ~ x, cut = cut)
  )
  expect_equal(unname(coef(fic$fits[[1]])), unname(coef(ref)), tolerance = 1e-6)
})

test_that("pooled add_* produce ordered, monotone quantities", {
  icd <- make_ic_data(250, seed = 2)
  cut <- seq(0, 10, by = 0.5)
  fit <- pamm_ic(Surv(L, R, type = "interval2") ~ x, icd, cut = cut, m = 6)
  ped <- as_ped(icd, Surv(L, R, type = "interval2") ~ x, cut = cut)
  nd <- make_newdata(ped, tend = unique(tend))

  s <- add_surv_prob(nd, fit, nsim = 200)
  cu <- add_cumu_hazard(nd, fit, nsim = 200)
  h <- add_hazard(nd, fit, nsim = 200)

  expect_true(all(c("surv_prob", "surv_lower", "surv_upper") %in% names(s)))
  expect_true(all(diff(s$surv_prob) <= 1e-8)) # monotone
  expect_true(all(s$surv_lower <= s$surv_prob + 1e-8))
  expect_true(all(s$surv_prob <= s$surv_upper + 1e-8))
  expect_true(all(diff(cu$cumu_hazard) >= -1e-8))
  expect_true(all(
    h$ci_lower <= h$hazard + 1e-8 & h$hazard <= h$ci_upper + 1e-8
  ))

  nd_sparse <- make_newdata(ped, tend = c(5, 10))
  s_sparse <- add_surv_prob(nd_sparse, fit, ci = FALSE)
  s_full <- add_surv_prob(nd, fit, ci = FALSE)
  expect_equal(
    s_sparse$surv_prob,
    s_full$surv_prob[match(s_sparse$tend, s_full$tend)],
    tolerance = 1e-12
  )

  set.seed(202)
  cu_sparse <- add_cumu_hazard(nd_sparse, fit, nsim = 80)
  set.seed(202)
  cu_full <- add_cumu_hazard(nd, fit, nsim = 80)
  full_ix <- match(cu_sparse$tend, cu_full$tend)
  expect_equal(
    cu_sparse$cumu_lower,
    cu_full$cumu_lower[full_ix],
    tolerance = 1e-12
  )
  expect_equal(
    cu_sparse$cumu_upper,
    cu_full$cumu_upper[full_ix],
    tolerance = 1e-12
  )
})

test_that("competing-risks IC pipeline yields valid pooled CIFs", {
  set.seed(5)
  df <- data.frame(x = runif(300, -1, 1))
  sdf <- sim_pexp(~ -2.5 + 0.6 * x, df, cut = seq(0, 10, by = 0.25))
  sdf$cause <- ifelse(
    sdf$status == 1,
    sample(c(1, 2), nrow(sdf), TRUE, c(.6, .4)),
    0
  )
  icd <- add_inspections(sdf, rate = 1, max_time = 10)
  cut <- seq(0, 10, by = 0.5)

  fcr <- pamm_ic_cr(
    Surv(L, R, type = "interval2") ~ x,
    icd,
    cause = "cause",
    cut = cut,
    m = 4
  )
  expect_identical(fcr$type, "cr")
  expect_setequal(fcr$cause_levels, c("1", "2"))

  ped_cr <- suppressWarnings(as_ped(
    transform(
      icd,
      time = pmin(true_time, 10),
      status = ifelse(true_time > 10, 0, cause)
    ),
    Surv(time, status) ~ x,
    cut = cut
  ))
  nd <- make_newdata(ped_cr, tend = unique(tend), cause = unique(cause))
  nd <- dplyr::group_by(nd, cause)
  cif <- add_cif(nd, fcr, nsim = 120)

  expect_true(all(c("cif", "cif_lower", "cif_upper") %in% names(cif)))
  expect_true(all(cif$cif >= 0 & cif$cif <= 1))
  expect_true(all(tapply(
    cif$cif,
    cif$cause,
    function(z) all(diff(z) >= -1e-8)
  )))
  expect_true(all(
    cif$cif_lower <= cif$cif + 1e-8 &
      cif$cif <= cif$cif_upper + 1e-8
  ))

  cif1 <- add_cif(nd, fcr, ci = FALSE, nsim = 20)
  cif2 <- add_cif(nd, fcr, ci = FALSE, nsim = 120)
  expect_equal(cif1$cif, cif2$cif, tolerance = 1e-12)
  expect_false(any(c("cif_lower", "cif_upper") %in% names(cif1)))
  expect_error(add_cif(dplyr::ungroup(nd), fcr, ci = FALSE), "group by cause")
})

test_that("competing-risks imputation samples exact unknown causes", {
  set.seed(10)
  df <- data.frame(
    id = 1:30,
    x = runif(30, -1, 1),
    L = c(rep(1, 5), rep(2, 5), rep(0, 10), rep(3, 10)),
    R = c(rep(1, 5), rep(2, 5), rep(2, 10), rep(Inf, 10)),
    cause = c(
      1,
      2,
      NA,
      1,
      2,
      rep(c(1, 2), length.out = 5),
      NA,
      rep(c(1, 2), length.out = 9),
      rep(0, 10)
    )
  )
  cut <- seq(0, 4, by = 0.25)
  formula <- Surv(L, R, type = "interval2") ~ .
  ic <- parse_ic_surv(formula, df, id = "id")

  cause_raw <- df$cause
  is_cens <- !is.na(cause_raw) & as.character(cause_raw) == "0"
  is_cens <- is_cens | as.character(ic$ic_kind) == "right"
  cause_known <- ifelse(is_cens, NA, as.character(cause_raw))
  cause_levels <- sort(unique(stats::na.omit(cause_known)))

  L <- ic$ic_L
  R <- ic$ic_R
  t_mid <- pmin(ifelse(ic$ic_kind == "left", R / 2, (L + R) / 2), max(cut))
  cause0 <- cause_known
  unknown <- which(!is_cens & is.na(cause0))
  cause0[unknown] <- sample(cause_levels, length(unknown), replace = TRUE)

  ped0 <- build_ic_ped_cr(
    ic,
    t_mid,
    cause0,
    is_cens,
    cut,
    formula,
    "id",
    cause_levels,
    0,
    cause_var = "cause"
  )
  fit0 <- pamm(ped_status ~ s(tend, by = cause, k = 5) + cause + x, data = ped0)
  cache <- ic_pred_cache(fit0, ic, cut, cause_levels = cause_levels)
  imp <- impute_ic_cr(
    fit0,
    ic,
    cut,
    beta = coef(fit0),
    cache = cache,
    cause_known = cause_known
  )

  exact_unknown <- which(
    as.character(ic$ic_kind) == "exact" &
      is.na(cause_known) &
      !is_cens
  )
  expect_false(anyNA(imp$cause[exact_unknown]))
  expect_true(all(imp$cause[exact_unknown] %in% cause_levels))
})

test_that("known-cause CR imputation does not keep rejected proposals", {
  cut <- c(0, 1)
  ic <- data.frame(
    id = 1L,
    ic_L = 0,
    ic_R = 1,
    ic_kind = factor(
      "interval",
      levels = c("exact", "right", "left", "interval")
    )
  )
  cache <- list(
    ii = int_info(cut),
    n_int = 1L,
    n_sub = 1L,
    X_list = list(
      a = matrix(-Inf, nrow = 1L, ncol = 1L),
      b = matrix(0, nrow = 1L, ncol = 1L)
    ),
    cause_levels = c("a", "b")
  )

  expect_error(
    impute_ic_cr(
      list(),
      ic,
      cut,
      beta = 1,
      cache = cache,
      cause_known = "a"
    ),
    "Failed to draw an interval-censored event time for known cause"
  )
})

test_that("add_inspections brackets the true event time", {
  set.seed(7)
  df <- data.frame(x = runif(100, -1, 1))
  sdf <- sim_pexp(~ -2 + 0.3 * x, df, cut = seq(0, 10, by = 0.25))
  icd <- add_inspections(sdf, rate = 2, max_time = 10)
  ev <- icd$status == 1 & is.finite(icd$R)
  expect_true(all(icd$true_time[ev] > icd$L[ev] - 1e-9))
  expect_true(all(icd$true_time[ev] <= icd$R[ev] + 1e-9))
  # fixed schedule mechanism
  icf <- add_inspections(
    sdf,
    mechanism = "fixed",
    schedule = seq(1, 10, by = 1)
  )
  expect_true(all(icf$L %in% c(0, seq(1, 10, by = 1))))
})

test_that("pamm_ic iter > 1 performs chained re-imputation", {
  set.seed(21)
  df <- data.frame(x = runif(120, -1, 1))
  sdf <- sim_pexp(~ -2 + 0.4 * x, df, cut = seq(0, 10, by = 0.5))
  icd <- add_inspections(sdf, rate = 0.5, max_time = 10)
  ic_form <- Surv(L, R, type = "interval2") ~ x

  # default (iter = 1) is unchanged draw-for-draw by the refactor
  set.seed(1)
  f1 <- pamm_ic(ic_form, icd, m = 2, cut = seq(0, 10, by = 0.5))
  set.seed(1)
  f1b <- pamm_ic(ic_form, icd, m = 2, iter = 1, cut = seq(0, 10, by = 0.5))
  expect_equal(coef(f1$fits[[1]]), coef(f1b$fits[[1]]))
  expect_equal(f1$iter, 1)

  # iter = 2 runs, records iter, and yields different chain fits
  set.seed(1)
  f2 <- pamm_ic(ic_form, icd, m = 2, iter = 2, cut = seq(0, 10, by = 0.5))
  expect_s3_class(f2, "pamm_ic")
  expect_length(f2$fits, 2)
  expect_equal(f2$iter, 2)
  expect_false(isTRUE(all.equal(coef(f1$fits[[1]]), coef(f2$fits[[1]]))))
  # pooled machinery still works downstream
  ped <- as_ped(icd, ic_form, cut = seq(0, 10, by = 0.5))
  nd <- make_newdata(ped, tend = c(5, 10))
  s2 <- add_surv_prob(nd, f2, nsim = 50)
  expect_true(all(s2$surv_lower <= s2$surv_prob + 1e-8))
  expect_true(all(s2$surv_prob <= s2$surv_upper + 1e-8))
})

test_that("pamm_ic_cr supports iter and records it", {
  set.seed(31)
  n <- 150
  df <- data.frame(id = seq_len(n), x = runif(n, -1, 1))
  t1 <- rexp(n, 0.15)
  t2 <- rexp(n, 0.10)
  obs <- pmin(t1, t2, 10)
  df$time <- obs
  df$status <- as.integer(obs < 10)
  df$cause <- ifelse(obs >= 10, 0L, ifelse(t1 < t2, 1L, 2L))
  icd <- add_inspections(df, rate = 0.8, max_time = 10)
  set.seed(2)
  f2 <- pamm_ic_cr(
    Surv(L, R, type = "interval2") ~ x,
    icd,
    cause = "cause",
    cut = seq(0, 10, by = 1),
    m = 2,
    iter = 2
  )
  expect_s3_class(f2, "pamm_ic")
  expect_equal(f2$type, "cr")
  expect_equal(f2$iter, 2)
  expect_length(f2$fits, 2)
  expect_identical(f2$unstable_chains, integer(0))
})

test_that("warn_unstable_chains flags degenerate imputation fits", {
  ok_fit <- list(coefficients = c(a = -1.5, b = 0.4), Vp = diag(c(0.04, 0.01)))
  bad_coef <- list(coefficients = c(a = 46, b = 0.4), Vp = diag(c(0.04, 0.01)))
  bad_se <- list(coefficients = c(a = -1.5, b = 0.4), Vp = diag(c(1e30, 0.01)))
  expect_silent(
    expect_identical(warn_unstable_chains(list(ok_fit, ok_fit)), integer(0))
  )
  expect_warning(
    out <- warn_unstable_chains(list(ok_fit, bad_coef, bad_se)),
    "unstable"
  )
  expect_identical(out, c(2L, 3L))
})

test_that("add_inspections right-censoring is non-informative", {
  set.seed(11)
  df <- data.frame(x = runif(400, -1, 1))
  # high baseline + horizon 4 < max event time: plenty of survivors (status 0)
  sdf <- sim_pexp(~ -1 + 0.3 * x, df, cut = seq(0, 4, by = 0.25))
  expect_gt(sum(sdf$status == 0), 0)

  # terminal_exam = TRUE (default): every event before the horizon has a
  # finite upper bound, so right-censored records never hide events in
  # (L, max_time] -- the regression guard for the informative-censoring bug
  icd <- add_inspections(sdf, rate = 0.3, max_time = 4)
  ev <- icd$status == 1
  expect_true(all(is.finite(icd$R[ev & icd$true_time <= 4])))
  expect_false(any(is.infinite(icd$R) & icd$L < 4))
  expect_true(all(icd$L[icd$status == 0] == icd$true_time[icd$status == 0]))

  # terminal_exam = FALSE: survivors are censored at their LAST INSPECTION,
  # never at their exact exit time (the other coherent convention)
  icn <- add_inspections(sdf, rate = 0.3, max_time = 4, terminal_exam = FALSE)
  sv <- icn$status == 0
  expect_true(all(is.infinite(icn$R[sv])))
  expect_true(all(icn$L[sv] <= icn$true_time[sv]))
  # with Exp(0.3) gaps, last inspections a.s. differ from the exit time
  expect_true(all(icn$L[sv] < icn$true_time[sv]))

  # NA status is rejected, not silently treated as an event
  sdf_na <- sdf
  sdf_na$status[1] <- NA
  expect_error(add_inspections(sdf_na, rate = 1), "missing")
})

test_that("pooled pamm_ic adders warn on undergrouped newdata", {
  set.seed(11)
  df <- data.frame(z = c(rep(0, 200), rep(1, 200)))
  sdf <- sim_pexp(~ -2 + 1.0 * z, df, cut = seq(0, 8, by = 0.25))
  icd <- add_inspections(sdf, rate = 1, max_time = 8)
  cut <- seq(0, 8, by = 0.5)
  fit <- pamm_ic(Surv(L, R, type = "interval2") ~ s(tend) + z, icd, cut = cut, m = 2)
  ped <- as_ped(icd, Surv(L, R, type = "interval2") ~ z, cut = cut)

  # multiple covariate groups but no group_by() -> distinct curves get pooled
  nd_multi <- make_newdata(ped, tend = unique(tend), z = c(0, 1))
  expect_warning(add_surv_prob(nd_multi, fit, ci = FALSE), "group_by")
  expect_warning(add_cumu_hazard(nd_multi, fit, ci = FALSE), "group_by")

  # group_by() silences the warning and yields correct per-group curves
  nd_g <- dplyr::group_by(nd_multi, z)
  expect_warning(s <- add_surv_prob(nd_g, fit, ci = TRUE, nsim = 40), regexp = NA)
  s <- dplyr::ungroup(s)
  expect_true(all(tapply(s$surv_prob, s$z, function(v) all(diff(v) <= 1e-8))))
  expect_true(all(s$surv_lower <= s$surv_prob + 1e-8 & s$surv_prob <= s$surv_upper + 1e-8))

  # a single-curve prediction grid must not warn
  nd1 <- make_newdata(ped, tend = unique(tend))
  expect_warning(add_surv_prob(nd1, fit, ci = FALSE), regexp = NA)
})
