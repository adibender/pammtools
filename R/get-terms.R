#' Enumerate plottable univariate smooth terms of a fitted model
#'
#' Internal helper. Given the model \code{data} and a fitted \code{\link[mgcv]{gam}}
#' object, returns a tibble with one row per smooth \emph{curve} to be drawn by
#' \code{\link{get_terms}} / \code{\link{gg_smooth}}. Only smooths that vary over
#' exactly one \emph{numeric} covariate are returned. This includes ordinary 1d
#' smooths (\code{s()}, 1d \code{ti()}), by-variable smooths
#' (\code{s(x, by = z)}) and factor-smooth interactions
#' (\code{s(x, fac, bs = "fs")}, \code{s(x, fac, bs = "sz")}). Tensor and
#' multivariate smooths (\code{te()}, \code{t2()}, 2d \code{ti()}, \code{s(x, z)})
#' as well as (correlated) random effects (\code{bs = "re"}, \code{bs = "mrf"})
#' are excluded -- use \code{\link{gg_tensor}} / \code{\link{gg_re}} for those.
#'
#' Smooths that are indexed by a factor (a factor \code{by}-variable or the factor
#' in an \code{fs}/\code{sz} interaction) are expanded into one row per factor
#' level, all sharing the same \code{facet} so that \code{\link{gg_smooth}} can
#' draw them in a single panel, distinguished by colour/fill.
#'
#' Returns \code{NULL} for models without a \code{$smooth} component (e.g.
#' \code{\link[survival]{coxph}}), in which case \code{\link{get_terms}} falls back
#' to label-based extraction.
#'
#' @param data A data frame containing the variables used to fit the model.
#' @param fit A fitted model object.
#' @return A tibble with columns \code{facet}, \code{level}, \code{var},
#' \code{col} and the list-column \code{settings}, or \code{NULL}.
#' @keywords internal
get_smooth_terms <- function(data, fit) {
  smooths <- fit[["smooth"]]
  # NULL only for fits without mgcv smooth metadata (-> non-gam fallback). An
  # mgcv fit whose smooths are all excluded (e.g. only te()/random effects) still
  # has a `$smooth` component and should yield an *empty* selection here, not the
  # non-gam fallback.
  if (is.null(smooths)) {
    return(NULL)
  }

  rows <- purrr::compact(purrr::map(smooths, smooth_term_rows, data = data))
  if (length(rows) == 0) {
    return(empty_smooth_terms())
  }
  dplyr::bind_rows(rows)
}

#' @keywords internal
empty_smooth_terms <- function() {
  tibble::tibble(
    facet    = character(),
    level    = character(),
    var      = character(),
    col      = character(),
    settings = list()
  )
}

#' Turn a single mgcv smooth into zero or more curve specifications
#'
#' @param s A single smooth object from \code{fit$smooth}.
#' @inheritParams get_smooth_terms
#' @keywords internal
smooth_term_rows <- function(s, data) {
  # (correlated) random effects are not "x vs f(x)" curves -> handled by gg_re()
  if (inherits(s, c("random.effect", "mrf.smooth"))) {
    return(NULL)
  }

  term_vars <- s[["term"]]
  is_num <- vapply(
    term_vars,
    function(v) v %in% names(data) && is.numeric(data[[v]]),
    logical(1)
  )
  num_vars <- term_vars[is_num]
  fac_vars <- term_vars[!is_num]

  # only smooths over exactly one numeric covariate are 1d-plottable.
  # NB: use the number of *numeric* covariates, not s$dim, because factor-smooth
  # interactions (bs = "sz") may report dim > 1 while varying over a single
  # numeric covariate.
  if (length(num_vars) != 1) {
    return(NULL)
  }
  var <- num_vars
  label <- s[["label"]]

  by <- s[["by"]]
  has_by <- !is.null(by) && !is.na(by) && by != "NA"
  by_is_factor <- has_by && by %in% names(data) && is_categorical(data[[by]])

  # Case A: factor by-variable. mgcv creates one smooth (and one term column) per
  # level, each carrying its own `by.level`. Collapse them into a shared facet.
  if (by_is_factor) {
    lvl <- s[["by.level"]]
    if (is.null(lvl)) {
      lvl <- parse_by_level(label, by, levels(as.factor(data[[by]])))
    }
    lvl <- as.character(lvl)
    facet <- substr(label, 1, nchar(label) - nchar(lvl))
    return(tibble::tibble(
      facet    = facet,
      level    = lvl,
      var      = var,
      col      = label,
      settings = list(stats::setNames(list(lvl), by))
    ))
  }

  # Case B: factor-smooth interaction (bs = "fs"/"sz"). A single term column whose
  # value depends on the level of the factor covariate set in `newdata`.
  if (length(fac_vars) >= 1) {
    fac  <- fac_vars[1]
    lvls <- levels(as.factor(data[[fac]]))
    return(tibble::tibble(
      facet    = label,
      level    = lvls,
      var      = var,
      col      = label,
      settings = lapply(lvls, function(l) stats::setNames(list(l), fac))
    ))
  }

  # Case C: plain numeric smooth, optionally with a numeric by-variable. For a
  # numeric by-variable the term contribution is f(var) * by, so we set by = 1 to
  # recover the base smooth f(var).
  setting <- list()
  if (has_by) {
    setting <- stats::setNames(list(1), by)
  }
  tibble::tibble(
    facet    = label,
    level    = NA_character_,
    var      = var,
    col      = label,
    settings = list(setting)
  )
}

#' @keywords internal
is_categorical <- function(x) is.factor(x) || is.character(x) || is.logical(x)

#' Parse the factor level from a factor-by smooth label
#'
#' Fallback for the rare case where \code{by.level} is unavailable. A label such as
#' \code{"s(tend):metastasesyes"} encodes the level (\code{"yes"}) as the suffix of
#' the by-variable name (\code{"metastases"}).
#' @keywords internal
parse_by_level <- function(label, by, lvls) {
  suffix <- sub("^.*:", "", label)
  lvl <- if (startsWith(suffix, by)) substring(suffix, nchar(by) + 1L) else suffix
  if (!lvl %in% lvls) {
    hit <- lvls[vapply(lvls, function(l) endsWith(label, l), logical(1))]
    if (length(hit) >= 1) lvl <- hit[which.max(nchar(hit))]
  }
  if (!lvl %in% lvls) {
    stop(
      sprintf(
        "Could not determine the factor level for smooth '%s' (by = '%s').",
        label, by
      ),
      call. = FALSE
    )
  }
  lvl
}

#' Resolve requested terms against the model's plottable smooths
#'
#' @param smooth_tbl Output of \code{\link{get_smooth_terms}}.
#' @param terms A character vector of requested terms, or \code{NULL} for all.
#' @keywords internal
resolve_terms <- function(smooth_tbl, terms) {
  if (is.null(terms)) {
    return(smooth_tbl)
  }

  matched <- purrr::map(terms, function(t) {
    # exact match on the smooth label / facet first, then the (numeric) variable
    hit <- smooth_tbl[smooth_tbl$facet == t | smooth_tbl$col == t, , drop = FALSE]
    if (nrow(hit) == 0) {
      hit <- smooth_tbl[smooth_tbl$var == t, , drop = FALSE]
    }
    if (nrow(hit) == 0) {
      warning(
        sprintf(
          "'%s' is not a univariate smooth term in the model and will be skipped.",
          t
        ),
        call. = FALSE
      )
    }
    hit
  })

  out <- dplyr::bind_rows(matched)
  out[!duplicated(out[c("facet", "level", "col")]), , drop = FALSE]
}

#' Extract the partial effect of a single smooth curve
#'
#' @param data A data frame containing variables used to fit the model. The first
#' row is used as the basis for all covariates other than the one being varied
#' (their values are irrelevant for the term-wise contribution).
#' @param fit A fitted object of class \code{\link[mgcv]{gam}}.
#' @param spec A single-row tibble (one row of \code{\link{get_smooth_terms}})
#' describing the curve to extract.
#' @param n Number of points at which to evaluate the smooth over the range of
#' its covariate.
#' @param conf_level The confidence level for the pointwise confidence interval.
#' @param ... Further arguments (currently unused).
#' @import dplyr
#' @importFrom stats predict
#' @keywords internal
get_term <- function(data, fit, spec, n = 100, conf_level = 0.95, ...) {
  z_value <- normal_ci_multiplier(conf_level)

  var      <- spec[["var"]]
  col      <- spec[["col"]]
  facet    <- spec[["facet"]]
  level    <- spec[["level"]]
  settings <- spec[["settings"]][[1]]

  # values at which the term contribution will be evaluated
  seq_term <- seq_range(pull(data, var), n = n)

  # first row as basis (values of other covariates irrelevant for type = "terms")
  new_df <- data[rep(1, length(seq_term)), , drop = FALSE]
  new_df[[var]] <- seq_term

  # activate the relevant factor level / by-variable value
  for (nm in names(settings)) {
    val <- settings[[nm]]
    if (is_categorical(data[[nm]])) {
      new_df[[nm]] <- factor(val, levels = levels(as.factor(data[[nm]])))
    } else {
      new_df[[nm]] <- as.numeric(val)
    }
  }

  term_info <- predict(fit, newdata = new_df, type = "terms", se.fit = TRUE)
  idx <- match(col, colnames(term_info$fit))
  if (is.na(idx)) {
    stop(
      sprintf("Could not find term column '%s' in model predictions.", col),
      call. = FALSE
    )
  }

  tibble::tibble(
    term  = facet,
    x     = seq_term,
    level = as.character(level),
    eff   = as.numeric(term_info$fit[, idx]),
    se    = as.numeric(term_info$se.fit[, idx])
  ) %>%
    mutate(
      ci_lower = .data$eff - z_value * .data$se,
      ci_upper = .data$eff + z_value * .data$se
    )
}

#' Extract a partial effect for models without a \code{$smooth} component
#'
#' Fallback used for fits such as \code{\link[survival]{coxph}} that support
#' \code{predict(type = "terms")} but expose no mgcv smooth metadata. Matching is
#' anchored to the variable name (exact, or as a parenthesised argument such as
#' \code{pspline(karno)}) rather than an unanchored substring.
#'
#' @inheritParams get_term
#' @param term A character string naming the model term/variable.
#' @keywords internal
get_term_legacy <- function(data, fit, term, n = 100, conf_level = 0.95, ...) {
  z_value <- normal_ci_multiplier(conf_level)

  seq_term <- seq_range(pull(data, term), n = n)
  new_df <- data[rep(1, length(seq_term)), , drop = FALSE]
  new_df[[term]] <- seq_term

  term_info <- predict(fit, newdata = new_df, type = "terms", se.fit = TRUE)
  cn <- colnames(term_info$fit)

  idx <- which(cn == term)
  if (length(idx) == 0) {
    # match `term` as a whole "word" (e.g. as the argument of pspline(karno)),
    # not as an arbitrary substring. `\Q...\E` literal quoting requires PCRE.
    idx <- grep(
      paste0("(^|[^[:alnum:]._])\\Q", term, "\\E($|[^[:alnum:]._])"),
      cn,
      perl = TRUE
    )
  }
  if (length(idx) == 0) {
    stop(
      paste0("No model term matched `term`: ", term),
      call. = FALSE
    )
  }
  idx <- idx[1]

  tibble::tibble(
    term  = term,
    x     = seq_term,
    level = NA_character_,
    eff   = as.numeric(term_info$fit[, idx]),
    se    = as.numeric(term_info$se.fit[, idx])
  ) %>%
    mutate(
      ci_lower = .data$eff - z_value * .data$se,
      ci_upper = .data$eff + z_value * .data$se
    )
}

#' @keywords internal
empty_term_tibble <- function() {
  tibble::tibble(
    term     = character(),
    x        = double(),
    level    = character(),
    eff      = double(),
    se       = double(),
    ci_lower = double(),
    ci_upper = double()
  )
}

#' Extract the partial effects of univariate smooth model terms
#'
#' Creates, for each requested univariate smooth, a sequence over the range of the
#' smooth's numeric covariate, evaluates the term-wise contribution via
#' \code{predict(fit, newdata = ., type = "terms")} and stacks the results into a
#' tidy data frame.
#'
#' For \code{\link[mgcv]{gam}} fits the requested \code{terms} are matched against
#' the model's smooths (see \code{\link{get_smooth_terms}}): a bare variable name
#' (e.g. \code{"tend"}) selects \emph{every} univariate smooth over that variable
#' -- the main effect \code{s(tend)} as well as any \code{s(tend, by = ...)} or
#' factor-smooth interaction -- while an exact smooth label (e.g. \code{"s(tend)"})
#' selects a single smooth. Names that do not match any smooth (for example
#' parametric factor main effects) are skipped with a warning; use
#' \code{\link{gg_fixed}} for those. For factor-indexed smooths one curve per
#' factor level is returned, identified by the \code{level} column.
#'
#' For models without mgcv smooth metadata (e.g. \code{\link[survival]{coxph}})
#' \code{terms} must be supplied and is matched against the columns of
#' \code{predict(type = "terms")}.
#'
#' @inheritParams get_term
#' @param terms A character vector (can be length one) specifying the terms for
#' which partial effects will be returned. If \code{NULL} (the default) all
#' univariate smooth terms in the model are used.
#' @param ... Further arguments controlling extraction, passed on per term, e.g.
#' \code{n} (number of evaluation points) and \code{conf_level}.
#' @import checkmate
#' @importFrom purrr map_dfr map compact
#' @return A tibble with columns \code{term}, \code{x}, \code{level}, \code{eff},
#' \code{se}, \code{ci_lower} and \code{ci_upper}.
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ pspline(karno) + pspline(age), data=veteran)
#' terms_df <- veteran %>% get_terms(fit, terms = c("karno", "age"))
#' head(terms_df)
#' tail(terms_df)
#' @export
get_terms <- function(data, fit, terms = NULL, ...) {
  # check inputs
  assert_class(data, "data.frame")
  if (!is.null(terms)) {
    assert_character(terms, min.len = 1, unique = TRUE)
  }

  smooth_tbl <- get_smooth_terms(data, fit)

  # fallback for models without mgcv smooth metadata (e.g. coxph)
  if (is.null(smooth_tbl)) {
    if (is.null(terms)) {
      stop(
        "`terms` must be supplied for models without a `$smooth` component.",
        call. = FALSE
      )
    }
    return(map_dfr(terms, function(t) get_term_legacy(data, fit, t, ...)))
  }

  selected <- resolve_terms(smooth_tbl, terms)
  if (nrow(selected) == 0) {
    warning(
      "No matching univariate smooth terms found; returning an empty result.",
      call. = FALSE
    )
    return(empty_term_tibble())
  }

  map_dfr(
    seq_len(nrow(selected)),
    function(i) get_term(data = data, fit = fit, spec = selected[i, ], ...)
  )
}
