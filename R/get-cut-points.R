#' Obtain interval break points
#'
#'
#' @import Formula
#' @keywords internal
get_cut <- function(data, formula, cut=NULL, ...) {
  UseMethod("get_cut", data)
}

#' @rdname get_cut
#' @inherit get_cut
get_cut.default <- function(data, formula, cut=NULL, max_time=NULL, ...) {

  outcom_vars <- get_lhs_vars(formula)
  if (is.null(cut)) {
    cut <- unique(data[["ped_time"]][data[["ped_status"]] == 1])
    if(!is.null(max_time)) {
      cut <- cut[cut < max_time]
      cut <- c(cut, max_time)
    }
  }

  # sort interval cut points in case they are not (so that interval factor
  # variables will be in correct ordering)
  sort(cut)

}

# get_cut.nested_fdf <- function(data, formula, cut=NULL, max_time=NULL, ...) {


# }
