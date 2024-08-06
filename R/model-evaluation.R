#' Transform crps object to data.frame
#'
#' A\code{as.data.frame} S3 method for objects of class \code{\link[pec]{crps}}.
#'
#' @inheritParams base::as.data.frame
#' @param x An object of class \code{crps}. See \code{\link[pec]{crps}}.
#' @importFrom tidyr pivot_longer
#'
#' @export
as.data.frame.crps <- function(x, row.names = NULL, optional = FALSE, ...) {

  m <- matrix(x, nrow = dim(x)[1], ncol = dim(x)[2])
  colnames(m) <- attr(x, "dimnames")[[2]]

  m <- as.data.frame(m)
  m$method <- attr(x, "dimnames")[[1]]

  m <- m %>%
    pivot_longer(cols = -one_of("method"), values_to = "IBS") %>%
    dplyr::rename(time = "name")

}
