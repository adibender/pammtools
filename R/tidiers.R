## Functions to extract effect information from fitted models

normal_ci_multiplier <- function(level = 0.95) {
  stats::qnorm(1 - (1 - level) / 2)
}


#' Calculate confidence intervals
#'
#' Given 2 column matrix or data frame, returns 3 column data.frame
#' with coefficient estimate plus lower and upper borders of the
#' 95% confidence intervals.
#'
#' @param ftab A table with two columns, containing coefficients in the first
#' column and standard-errors in the second column.
#' @importFrom tibble as_tibble
#' @keywords internal
calc_ci <- function(ftab) {
  colnames(ftab) <- c("coef", "se")
  rnames <- rownames(ftab)
  ci_mult <- normal_ci_multiplier()
  ftab <- as_tibble(ftab)
  ftab$variable <- rnames
  ftab$ci_lower <- ftab$coef - ci_mult * ftab$se
  ftab$ci_upper <- ftab$coef + ci_mult * ftab$se
  ftab$se <- NULL

  ftab[, c("variable", "coef", "ci_lower", "ci_upper")]
}

#' Extract fixed coefficient table from model object
#'
#' Given a model object, returns a data frame with columns \code{variable},
#' \code{coef} (coefficient), \code{ci_lower} (lower 95\% CI) and
#' \code{ci_upper} (upper 95\% CI).
#'
#' @param x A model object.
#' @param ... Currently not used.
#' @export
tidy_fixed <- function(x, ...) {
  UseMethod("tidy_fixed", x)
}

#' @rdname tidy_fixed
#' @param intercept Should intercept also be returned? Defaults to \code{FALSE}.
#' @export
tidy_fixed.gam <- function(x, intercept = FALSE, ...) {
  ftab <- summary(x)[["p.table"]][, 1:2]
  if (!intercept) {
    ftab <- ftab[!grepl("Intercept", rownames(ftab)), , drop = FALSE]
  }
  calc_ci(ftab)
}

#' @rdname tidy_fixed
#' @export
tidy_fixed.scam <- function(x, intercept = FALSE, ...) {
  tidy_fixed.gam(x, intercept = intercept, ...)
}

#' @rdname tidy_fixed
#' @importFrom tibble as_tibble
#' @keywords internal
#' @examples
#' library(survival)
#' gc <- coxph(Surv(days, status)~age + sex, data = tumor)
#' tidy_fixed(gc)
#' @export
tidy_fixed.coxph <- function(x, ...) {
  ftab <- summary(x)[["coefficients"]][, c(1, 3)]
  calc_ci(ftab)
}


#' Extract 1d smooth objects in tidy data format.
#'
#' @rdname tidiers
#' @inheritParams get_plotinfo
#' @param keep A vector of variables to keep.
#' @param ci A logical value indicating whether confidence intervals should be
#' calculated and returned. Defaults to \code{TRUE}.
#' @param conf_level Numeric scalar in (0, 1). Confidence level used for the
#' returned confidence intervals when \code{ci = TRUE}. Defaults to
#' \code{0.95}.
#' @importFrom dplyr bind_rows
#' @export
tidy_smooth <- function(
  x,
  keep = c("x", "fit", "se", "xlab", "ylab"),
  ci = TRUE,
  conf_level = 0.95,
  ...
) {
  checkmate::assert_number(conf_level, lower = 0, upper = 1)
  if (conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be strictly between 0 and 1.", call. = FALSE)
  }
  z_value <- normal_ci_multiplier(conf_level)
  po <- get_plotinfo(x, ...)
  # index of list elements that are 1d smooths and not random effects
  ind1d <- vapply(
    X = po,
    FUN = function(z) !is.null(z[["x"]]) & is.null(z[["main"]]),
    FUN.VALUE = logical(1)
  )
  # keep only variables of interest
  po <- lapply(po[ind1d], "[", i = keep, drop = TRUE)

  # transform to data.frame
  po <- lapply(po, function(sm) {
    sm[["fit"]] <- as.vector(sm[["fit"]])
    temp <- as_tibble(sm)
    if (ci) {
      temp <- temp %>%
        mutate(
          ci_lower = .data$fit - z_value * .data$se,
          ci_upper = .data$fit + z_value * .data$se
        )
    }
    temp
  })

  return(bind_rows(po))
}


#' Extract 2d smooth objects in tidy format.
#'
#' @inheritParams tidy_smooth
#' @importFrom purrr cross_df
#' @importFrom tibble as_tibble
#' @import dplyr
#' @export
tidy_smooth2d <- function(
  x,
  keep = c("x", "y", "fit", "se", "xlab", "ylab", "main"),
  ci = FALSE,
  conf_level = 0.95,
  ...
) {
  checkmate::assert_number(conf_level, lower = 0, upper = 1)
  if (conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be strictly between 0 and 1.", call. = FALSE)
  }
  z_value <- normal_ci_multiplier(conf_level)
  po <- get_plotinfo(x, ...)

  ind2d <- vapply(
    X = po,
    FUN = function(z) !is.null(z[["x"]]) & !is.null(z[["y"]]),
    FUN.VALUE = logical(1)
  )

  # keep only variables of interest
  po <- lapply(po[ind2d], "[", i = keep, drop = TRUE)

  # transform to data.frame
  po <- lapply(po, function(sm) {
    sm[["fit"]] <- as.vector(sm[["fit"]])
    p1 <- as_tibble(sm[setdiff(keep, c("x", "y"))])
    xy <- cross_df(sm[c("x", "y")])
    xy <- bind_cols(xy, p1)
    if (ci) {
      xy <- xy %>%
        mutate(
          ci_lower = .data$fit - z_value * .data$se,
          ci_upper = .data$fit + z_value * .data$se
        )
    }
    xy
  })

  return(bind_rows(po))
}


#' Extract random effects in tidy data format.
#'
#' @inheritParams tidy_smooth
#' @importFrom dplyr bind_rows
#' @importFrom stats ppoints qnorm quantile
#' @rdname tidy_smooth
#' @seealso \code{\link[stats]{qqline}}
#' @export
tidy_re <- function(x, keep = c("fit", "main", "xlab", "ylab"), ...) {
  po <- get_plotinfo(x, ...)
  ind.re <- vapply(
    X = po,
    FUN = function(z) {
      (!is.null(z[["main"]])) & (z[["xlab"]] == "Gaussian quantiles")
    },
    FUN.VALUE = logical(1)
  )

  po <- lapply(po[ind.re], "[", i = keep, drop = TRUE)
  po <- lapply(po, function(z) {
    re.df <- do.call(cbind.data.frame, c(z, stringsAsFactors = FALSE))
    re.df$x <- qnorm(ppoints(length(re.df$fit))[order(order(re.df$fit))])
    # code to calculate qqslope and qqintercept from ?stats::qqline
    yl <- quantile(re.df$fit, probs = c(0.25, 0.75), type = 7, names = FALSE)
    xl <- qnorm(c(0.25, 0.75))
    re.df$qqslope <- diff(yl) / diff(xl)
    re.df$qqintercept <- yl[1L] - re.df$qqslope * xl[1L]

    re.df
  })

  return(bind_rows(po))
}


#' Extract plot information for all special model terms
#'
#' Given a \code{mgcv} \code{\link[mgcv]{gamObject}} (or a
#' \code{\link[scam]{scam}} object), returns the information
#' used for the default plots produced by \code{\link[mgcv]{plot.gam}}
#' (\code{\link[scam]{plot.scam}}, respectively).
#'
#' @inheritParams mgcv::plot.gam
#' @param ... Further arguments passed to \code{\link[mgcv]{plot.gam}}
#' @import mgcv
#' @importFrom checkmate assert_multi_class
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot
#' @export
get_plotinfo <- function(x, ...) {
  assert_multi_class(x, c("gam", "scam"))

  tmp <- paste0(tempfile(), ".png")
  png(tmp)
  po <- plot(x, page = 1, ...)
  dev.off()
  if (file.exists(tmp)) file.remove(tmp)

  class(po) <- c("mgcv.plotlist", class(po))

  return(po)
}
