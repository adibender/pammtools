#' Slim down a fitted PAMM for storage inside a \code{pamm_ic} object
#'
#' Removes the per-observation slots (model frame, fitted values, residuals,
#' working weights, ...) and the \code{call} (which captures the full PED data),
#' none of which are needed for the downstream multiple-imputation pooling: the
#' pooled \code{add_*} methods only require each fit's \code{coefficients},
#' \code{Vp}/\code{Ve} and the smooth/parametric structure used by
#' \code{predict(type = "lpmatrix")}. Stripping makes the stored size independent
#' of the data set size, so memory does not blow up with many imputations.
#'
#' @param fit A fitted \code{pamm}/\code{gam} object.
#' @return The same object with large per-observation slots removed; class and
#'   everything needed for \code{predict}/\code{coef}/\code{vcov} are retained.
#' @keywords internal
strip_pamm_fit <- function(fit) {
  drop <- c(
    "model",
    "y",
    "residuals",
    "fitted.values",
    "linear.predictors",
    "weights",
    "prior.weights",
    "offset",
    "wt",
    "working.weights",
    "hat",
    "z",
    "w",
    "std.rsd",
    "na.action"
  )
  for (s in intersect(drop, names(fit))) fit[[s]] <- NULL
  fit[["call"]] <- NULL
  fit
}

# Fit one imputation, capture its summary (for median-p pooling) while the full
# object is available, then return a stripped fit + PED row count. The full
# (unstripped) fit is also returned so the caller can use one model frame as a
# common training-grid reference for term-level MI diagnostics; only the slimmed
# fits are stored per imputation.
fit_strip_summarise <- function(model_formula, data, engine, ...) {
  f <- pamm(model_formula, data = data, engine = engine, ...)
  list(summary = summary(f), fit = strip_pamm_fit(f), full = f, n = nrow(data))
}

# Median over imputations of a per-fit table column (the "median-p rule" for
# pooling significance across MI, shown to work well for GAMs by Bolt et al.
# 2022, BMC Med Res Methodol, doi:10.1186/s12874-022-01613-w). Other columns
# (edf, statistics) are averaged.
pool_param_table <- function(ptabs, est, V, nsdf) {
  if (nsdf < 1L || is.null(ptabs[[1]]) || nrow(ptabs[[1]]) == 0L) return(NULL)
  rn <- rownames(ptabs[[1]])
  pcol <- ncol(ptabs[[1]])
  pmat <- vapply(ptabs, function(t) t[, pcol], numeric(nrow(ptabs[[1]])))
  pmed <- if (is.matrix(pmat)) apply(pmat, 1, stats::median) else
    stats::median(pmat)
  ind <- seq_len(nsdf)
  Est <- est[ind]
  SE <- sqrt(diag(V))[ind]
  out <- cbind(
    "Estimate" = Est,
    "Std. Error" = SE,
    "z value" = Est / SE,
    "Pr(>|z|)" = pmed
  )
  rownames(out) <- rn
  out
}

pool_smooth_table <- function(stabs) {
  if (is.null(stabs[[1]]) || nrow(stabs[[1]]) == 0L) return(NULL)
  rn <- rownames(stabs[[1]])
  arr <- simplify2array(stabs) # rows x cols x m
  if (length(dim(arr)) == 2L) arr <- array(arr, c(dim(arr), 1L))
  pc <- dim(arr)[2]
  out <- cbind(
    "edf" = apply(arr[, 1, , drop = FALSE], 1, mean),
    "Ref.df" = apply(arr[, 2, , drop = FALSE], 1, mean),
    "statistic" = apply(arr[, 3, , drop = FALSE], 1, mean),
    "p-value" = apply(arr[, pc, , drop = FALSE], 1, stats::median)
  )
  rownames(out) <- rn
  out
}

mi_riv <- function(within, between, m) {
  out <- within
  out[] <- NA_real_
  if (m <= 1L) {
    return(out)
  }

  between <- pmax(between, 0)
  ok <- is.finite(within) & within > 0 & is.finite(between)
  out[ok] <- (1 + 1 / m) * between[ok] / within[ok]

  no_var <- is.finite(within) & within <= 0 & between == 0
  out[no_var] <- 0
  only_between <- is.finite(within) & within <= 0 & between > 0
  out[only_between] <- Inf
  out
}

mi_fmi <- function(riv, m) {
  out <- riv
  out[] <- NA_real_
  if (m <= 1L) {
    return(out)
  }

  zero <- is.finite(riv) & riv == 0
  out[zero] <- 0

  finite <- is.finite(riv) & riv > 0
  dfbr <- (m - 1) * (1 + 1 / riv[finite])^2
  out[finite] <- (riv[finite] + 2 / (dfbr + 3)) / (riv[finite] + 1)

  out[is.infinite(riv) & riv > 0] <- 1
  out
}

pool_param_fmi_table <- function(fmi, nsdf) {
  if (nsdf < 1L || !length(fmi)) return(NULL)
  ind <- seq_len(nsdf)
  if (!any(is.finite(fmi[ind]))) return(NULL)
  out <- matrix(fmi[ind], ncol = 1L)
  colnames(out) <- "FMI"
  rownames(out) <- names(fmi)[ind]
  out
}

fmi_fivenum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(rep(NA_real_, 5L))
  }
  stats::quantile(
    x,
    probs = c(0, 0.25, 0.5, 0.75, 1),
    names = FALSE,
    na.rm = TRUE
  )
}

smooth_term_qoi <- function(fit, refdata, smooth_names) {
  labels <- vapply(fit[["smooth"]], `[[`, character(1L), "label")
  smooth_ix <- match(smooth_names, labels)
  if (anyNA(smooth_ix)) {
    stop("Smooth terms differ across imputation fits.", call. = FALSE)
  }

  X <- predict.gam(fit, refdata, type = "lpmatrix")
  q <- u <- matrix(
    NA_real_,
    nrow = nrow(X),
    ncol = length(smooth_names),
    dimnames = list(NULL, smooth_names)
  )
  beta <- coef(fit)

  for (kk in seq_along(smooth_names)) {
    sm <- fit[["smooth"]][[smooth_ix[kk]]]
    ind <- sm[["first.para"]]:sm[["last.para"]]
    Xs <- X[, ind, drop = FALSE]
    Xs <- sweep(Xs, 2L, colMeans(Xs), `-`)
    V <- fit[["Vp"]][ind, ind, drop = FALSE]

    q[, kk] <- as.numeric(Xs %*% beta[ind])
    u[, kk] <- pmax(rowSums((Xs %*% V) * Xs), 0)
  }

  list(q = q, u = u)
}

pool_smooth_fmi_table <- function(fits, skeleton, smooth_names) {
  if (
    length(fits) <= 1L ||
      is.null(skeleton) ||
      is.null(skeleton[["model"]]) ||
      !length(smooth_names)
  ) {
    return(NULL)
  }

  refdata <- skeleton[["model"]]
  mean_q <- m2_q <- W <- NULL
  used_terms <- NULL

  for (jj in seq_along(fits)) {
    qoi <- smooth_term_qoi(fits[[jj]], refdata, smooth_names)
    q <- qoi[["q"]]
    u <- qoi[["u"]]

    if (is.null(used_terms)) {
      used_terms <- colnames(q)
      mean_q <- m2_q <- W <- matrix(
        0,
        nrow = nrow(q),
        ncol = length(used_terms),
        dimnames = list(NULL, used_terms)
      )
    }
    if (!all(used_terms %in% colnames(q))) {
      stop("Smooth terms differ across imputation fits.", call. = FALSE)
    }

    q <- q[, used_terms, drop = FALSE]
    u <- u[, used_terms, drop = FALSE]

    W <- W + u
    delta <- q - mean_q
    mean_q <- mean_q + delta / jj
    m2_q <- m2_q + delta * (q - mean_q)
  }

  W <- W / length(fits)
  B <- m2_q / (length(fits) - 1)
  fmi <- mi_fmi(mi_riv(W, B, length(fits)), length(fits))
  out <- t(vapply(
    seq_len(ncol(fmi)),
    function(jj) {
      fmi_fivenum(fmi[, jj])
    },
    numeric(5L)
  ))
  colnames(out) <- c("Min", "Q1", "Median", "Q3", "Max")
  rownames(out) <- used_terms
  if (!any(is.finite(out))) {
    return(NULL)
  }
  out
}

#' Pool a list of (stripped) imputation fits into a pooled summary object
#'
#' Combines the \code{m} imputation fits with Rubin's rules: the pooled
#' parametric coefficients use \eqn{\bar Q} and
#' \eqn{V = \bar W + (1 + 1/m) B}. Smooth-term p-values are pooled with the
#' median-p rule (see references in \code{\link{strip_pamm_fit}}). Because
#' \code{mgcv} smooth coefficients can use different centered bases across
#' imputations, this returns a plain pooled summary object rather than a
#' \code{gam}: \code{add_*()} methods evaluate every imputation fit with its own
#' design matrix for predictions.
#'
#' @param fits List of stripped imputation fits.
#' @param smry List of \code{summary.gam} objects, one per fit (computed before
#'   stripping).
#' @param skeleton Optional full (unstripped) fit used only to supply a common
#'   training-grid model frame for smooth-term FMI summaries.
#' @keywords internal
#' @importFrom stats coef cov median
pool_pamm_fits <- function(fits, smry, skeleton = NULL) {
  m <- length(fits)
  p <- length(coef(fits[[1]]))
  cf <- vapply(fits, coef, numeric(p)) # p x m
  Qbar <- rowMeans(cf)
  nsdf <- fits[[1]][["nsdf"]]

  Wp <- Reduce(`+`, lapply(fits, function(f) f[["Vp"]])) / m
  We <- Reduce(
    `+`,
    lapply(fits, function(f) if (!is.null(f[["Ve"]])) f[["Ve"]] else f[["Vp"]])
  ) /
    m
  B <- if (m > 1) stats::cov(t(cf)) else matrix(0, p, p)
  Vp <- Wp + (1 + 1 / m) * B
  Ve <- We + (1 + 1 / m) * B

  riv <- mi_riv(diag(Wp), diag(B), m)
  fmi <- mi_fmi(riv, m)
  names(riv) <- names(fmi) <- names(Qbar)

  p.table <- pool_param_table(
    lapply(smry, `[[`, "p.table"),
    Qbar,
    Vp,
    nsdf
  )
  s.table <- pool_smooth_table(lapply(smry, `[[`, "s.table"))
  fmi.table <- pool_param_fmi_table(fmi, nsdf)
  smooth.fmi <- pool_smooth_fmi_table(
    fits,
    skeleton,
    rownames(s.table)
  )

  structure(
    list(
      family = fits[[1]][["family"]],
      nsdf = nsdf,
      p.table = p.table,
      s.table = s.table,
      riv = riv[seq_len(nsdf)],
      fmi = fmi[seq_len(nsdf)],
      fmi.table = fmi.table,
      smooth.fmi = smooth.fmi,
      Vp = Vp[seq_len(nsdf), seq_len(nsdf), drop = FALSE],
      Ve = Ve[seq_len(nsdf), seq_len(nsdf), drop = FALSE]
    ),
    class = "pamm_ic_pool"
  )
}

#' @rdname pamm_ic
#' @param x,object A \code{pamm_ic} object.
#' @param ... Passed on (ignored for \code{print}).
#' @export
print.pamm_ic <- function(x, ...) {
  cat("Piecewise exponential additive model for interval-censored data\n")
  cat(
    "  via multiple imputation (",
    x[["m"]],
    " proper imputations)\n",
    sep = ""
  )
  cat(
    "  task        : ",
    if (identical(x[["type"]], "cr")) "competing risks" else "single event",
    "\n",
    sep = ""
  )
  cat("  model       : ", deparse(x[["model_formula"]]), "\n", sep = "")
  n_ic <- sum(as.character(x[["ic"]][["ic_kind"]]) %in% c("interval", "left"))
  cat(
    "  subjects    : ",
    nrow(x[["ic"]]),
    " (",
    n_ic,
    " interval/left-censored)\n",
    sep = ""
  )
  cat(
    "  cut-points  : ",
    length(x[["cut"]]),
    " breaks in [",
    format(min(x[["cut"]])),
    ", ",
    format(max(x[["cut"]])),
    "]\n\n",
    sep = ""
  )
  cat(
    "Use summary() for pooled inference and add_*() for pooled",
    "quantities of interest.\n"
  )
  invisible(x)
}

#' @rdname pamm_ic
#' @export
summary.pamm_ic <- function(object, ...) {
  p <- object[["pooled"]]
  structure(
    list(
      type = object[["type"]],
      m = object[["m"]],
      model_formula = object[["model_formula"]],
      family = p[["family"]][["family"]],
      n_obs = object[["n_obs"]],
      n_subj = nrow(object[["ic"]]),
      n_ic = sum(
        as.character(object[["ic"]][["ic_kind"]]) %in%
          c("interval", "left")
      ),
      cut = object[["cut"]],
      p.table = p[["p.table"]],
      s.table = p[["s.table"]],
      fmi.table = p[["fmi.table"]],
      smooth.fmi = p[["smooth.fmi"]]
    ),
    class = "summary.pamm_ic"
  )
}

#' @rdname pamm_ic
#' @export
print.summary.pamm_ic <- function(x, ...) {
  cat(
    "Pooled PAMM summary (multiple imputation for interval-censored data)\n\n"
  )
  cat(
    "Task           :",
    if (identical(x[["type"]], "cr")) "competing risks" else "single event",
    "\n"
  )
  cat("Family         :", x[["family"]], "\n")
  cat("Model          :", deparse(x[["model_formula"]]), "\n")
  cat("Imputations    :", x[["m"]], "(proper)", "\n")
  cat(
    "Subjects       :",
    x[["n_subj"]],
    "(",
    x[["n_ic"]],
    "interval/left-censored ); PED rows per fit:",
    x[["n_obs"]],
    "\n\n"
  )

  if (!is.null(x[["p.table"]])) {
    cat("Parametric coefficients (Rubin-pooled estimates & SEs, median-p):\n")
    stats::printCoefmat(
      x[["p.table"]],
      has.Pvalue = TRUE,
      signif.stars = TRUE,
      digits = 3L
    )
  }
  if (!is.null(x[["s.table"]])) {
    cat(
      "\nApproximate significance of smooth terms",
      "(mean edf, median-p over imputations):\n"
    )
    stats::printCoefmat(
      x[["s.table"]],
      has.Pvalue = TRUE,
      signif.stars = TRUE,
      digits = 3L
    )
  }
  if (!is.null(x[["fmi.table"]]) || !is.null(x[["smooth.fmi"]])) {
    cat("\nFraction of missing information (FMI):\n")
  }
  if (!is.null(x[["fmi.table"]])) {
    cat("Parametric coefficients:\n")
    print(round(x[["fmi.table"]], 3L))
  }
  if (!is.null(x[["smooth.fmi"]])) {
    cat("\nSmooth terms (five-number summaries over training PED rows):\n")
    print(round(x[["smooth.fmi"]], 3L))
  }
  cat(
    "\nStandard errors include within- + between-imputation variance",
    "(Rubin's rules);\np-values are medians over imputations.\n"
  )
  invisible(x)
}
