#' Extract partial effects for specified model terms
#'
#' @param data A data frame containing variables used to fit the model. Only
#' first row will be used.
#' @param fit A fitted object of class \code{\link[mgcv]{gam}}.
#' @param term The (non-linear) model term of interest.
#' @param ... Further arguments passed to \code{\link{seq_range}}.
#' @inheritParams seq_range
#' @import dplyr
#' @importFrom stats predict
#' @importFrom rlang UQ
#' @keywords internal
get_term <- function(data, fit, term, n = 100, ...) {

  # values at which term contribution will be evaluated
  seq_term <- data %>% pull(term) %>% seq_range(n = n)

  # use first row as basis (values of other covariates irrelevant anyway)
  new_df <- data[1, ]

  # clean up as rest of the data not needed any longer
  rm(data)
  gc()

  term_name <- term
  # extract term contribution information (+ standard errors)
  new_df              <- new_df[rep(1, length(seq_term)), ]
  new_df[[term_name]] <- seq_term
  term_info           <- predict(fit, newdata = new_df, type = "terms",
    se.fit = TRUE)
  index_term          <- grep(term, colnames(term_info$fit), value = TRUE)

  new_df %>%
    mutate(
      term = term_name,
      eff  = as.numeric(term_info$fit[, index_term]),
      se   = as.numeric(term_info$se.fit[, index_term])) %>%
    mutate(
      ci_lower = .data$eff - 2 * .data$se,
      ci_upper = .data$eff + 2 * .data$se) %>%
  select(one_of(c("term", term_name, "eff", "se", "ci_lower", "ci_upper"))) %>%
  rename(x = UQ(term_name)) %>%
  as_tibble()

}

#' Extract the partial effects of non-linear model terms
#'
#' This function basically creates a new \code{df} from \code{data} for
#' each term in \code{terms}, creating a range from minimum and maximum of the
# respective terms. For each \code{df} it then calls
#' \code{predict(fit, newdata=df, type="terms")}. Terms are then
#' stacked to a tidy data frame.
#'
#' @inheritParams get_term
#' @param terms A character vector (can be length one). Specifies the terms
#' for which partial effects will be returned
#' @import checkmate
#' @importFrom purrr map_dfr
#' @return A tibble with 5 columns.
#' @examples
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ pspline(karno) + pspline(age), data=veteran)
#' terms_df <- veteran %>% get_terms(fit, terms = c("karno", "age"))
#' head(terms_df)
#' tail(terms_df)
#' @export
get_terms <- function(data, fit, terms, ...) {

  # check inputs
  assert_class(data, "data.frame")
  assert_character(terms, min.len = 1, unique = TRUE)

  # apply get_term to each element of terms
  map_dfr(terms, function(x) get_term(data = data, fit = fit, term = x), ...)

}
