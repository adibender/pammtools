.robust_fixture_cache <- new.env(parent = emptyenv())

get_tumor_ped_fixture <- function(n_rows = 120L) {
  key <- paste0("tumor_ped_", n_rows)
  if (!exists(key, envir = .robust_fixture_cache, inherits = FALSE)) {
    data("tumor", package = "pammtools")
    ped <- as_ped(
      data = dplyr::slice(tumor, seq_len(n_rows)),
      formula = survival::Surv(days, status) ~ age + complications,
      cut = c(0, 50, 100, 200, 300, 400),
      id = "id"
    )
    assign(key, ped, envir = .robust_fixture_cache)
  }

  get(key, envir = .robust_fixture_cache, inherits = FALSE)
}

get_multistate_fixture <- function(n_rows = 200L) {
  key <- paste0("multistate_", n_rows)
  if (!exists(key, envir = .robust_fixture_cache, inherits = FALSE)) {
    data("prothr", package = "mstate", envir = environment())
    prothr_subset <- prothr[seq_len(n_rows), ] |>
      dplyr::filter(.data$Tstart != .data$Tstop) |>
      dplyr::mutate(
        transition = as.factor(
          paste0(.data$from, "->", .data$to)
        )
      ) |>
      dplyr::select(-dplyr::one_of(c("trans", "treat")))

    ped_msm <- as_ped(
      data = prothr_subset,
      formula = survival::Surv(Tstart, Tstop, status) ~ .,
      transition = "transition",
      id = "id",
      timescale = "calendar"
    )

    pam_msm <- mgcv::gam(
      ped_status ~ s(tend, by = transition, bs = "cr") + transition,
      data = ped_msm,
      family = poisson(),
      offset = offset
    )

    assign(
      key,
      list(
        raw = prothr_subset,
        ped = ped_msm,
        pam = pam_msm
      ),
      envir = .robust_fixture_cache
    )
  }

  get(key, envir = .robust_fixture_cache, inherits = FALSE)
}

make_transition_newdata <- function(ped, shuffle = FALSE, seed = 42L) {
  ndf <- make_newdata(
    ped,
    tend = unique(tend),
    transition = unique(transition)
  ) |>
    dplyr::group_by(.data$transition) |>
    dplyr::arrange(.data$transition, .data$tend)

  if (shuffle) {
    set.seed(seed)
    ndf <- dplyr::slice_sample(ndf, prop = 1)
  }

  ndf
}

construct_fdl_smooth <- function(
  xt = list(),
  k = 8L,
  m = 2L,
  n = 50L,
  seed = 1L
) {
  set.seed(seed)
  design_data <- data.frame(x = runif(n))
  spec <- mgcv::s(
    x,
    bs = "ps",
    k = k,
    m = m,
    xt = xt
  )
  class(spec) <- "fdl.smooth.spec"

  pammtools:::smooth.construct.fdl.smooth.spec(
    object = spec,
    data = design_data,
    knots = NULL
  )
}
