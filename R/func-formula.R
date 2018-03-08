#' A formula special for defining functional covariates
#'
#' @rdname func
#' @importFrom purrr map
#' @export
#' @keywords internal
func <- function(..., by = NULL, ll_fun = function(t,te){te <= t}) {

  ## extract "variables" specified in ...
  vars   <- as.list(substitute(list(...)))[-1] %>% map(~as.character(.))
  by_var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  ## check if latency was specified (by checking if "-" contained in of of the ...)
  lgl_latency <- has_latency(vars)
  lgl_te      <- has_te(vars)
  lgl_time    <- has_time(vars)
  by_var      <- get_by_var(by_var)
  fun_covar   <- get_fun_covar(vars)
  if (is.null(by_var) & is.null(fun_covar)) {
    stop("No covariates specified.")
  }
  ## obtain name of functional covariate

  # return
  return(list(
    lgl_latency = lgl_latency,
    lgl_te      = lgl_te,
    lgl_time    = lgl_time,
    vars        = vars,
    by_var      = by_var,
    fun_covar   = fun_covar,
    ll_fun      = ll_fun))
}

#' @rdname func
#' @importFrom purrr map_lgl
#' @keywords internal
has_latency <- function(vars) {
    map_lgl(.x=vars, ~any(as.character(.x) == "-"))
}

#' @importFrom purrr map_lgl
#' @keywords internal
has_te <- function(vars) {
  map_lgl(.x=vars, ~any(as.character(.x)=="te")) & !has_latency(vars)
}

#' @importFrom purrr map_lgl
#' @keywords internal
has_time <- function(vars) {
  map_lgl(.x=vars, ~any(as.character(.x)=="t")) & !has_latency(vars)
}

#' @importFrom purrr invoke_map
#' @keywords internal
has_components <- function(vars) {
  lgl_components <- invoke_map(list(has_latency, has_te, has_time), vars=vars)
  names(lgl_components) <- c("lgl_latency", "lgl_te", "lgl_time")
  lgl_components
}

#' @inherit has_latency
#' @param vars Components exracted from \code{\link{func}} call.
#' @param exclude A character vector of components that should be viewed as
#' special components rather than variables in the data.
#' @importFrom purrr map_lgl flatten_chr
#' @keywords internal
get_fun_covar <- function(vars, exclude=c("t", "te", "-", "*")) {
  contains_time <- map_lgl(vars, ~any(. %in% c("t", "te")))
  fun_covar <- vars[!contains_time] %>% flatten_chr()
  # if(length(fun_covar > 1)) {
  #   stop("Only one, covariate (aside 't' and 'te' can be provided to 'func()'" )
  # }
  if(length(fun_covar) == 0) {
    fun_covar <- NULL
  }
  fun_covar
}

get_by_var <- function(by.var) {

  if (by.var == ".") stop("by=. not allowed")
  if (by.var == "NULL") {
    NULL
  } else{
    by.var
  }
}
