#' New basis for penalized lag selection
#'
#' Originally proposed in  Obermeier et al., 2015, Flexible Distributed Lags for Modelling Earthquake Data,
#' Journal of the Royal Statistical Society: Series C (Applied Statistics),
#' 10.1111/rssc.12077.
#' Here extended in order to penalize lead times in addition to lag times.
#' Ideally the lag-lead window would then be selected in a data-driven fashion.
#' Treat as experimental.
#'
#' @param object An object handled by mgcv
#' @param data The data set
#' @param knots A vector of knots
#' @keywords internal
smooth.construct.fdl.smooth.spec <- function(object, data, knots) {

  # modify object so that it's fitted as a p-spline signal regression term:
  object$bs <- "ps"
  object <- mgcv::smooth.construct.ps.smooth.spec(object, data, knots)

  if (!is.null(object$xt$fullrankpen) && object$xt$fullrankpen) {
    # add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions with
    # same variance as difference penalty : penalty = lambda * coef' (DiffPen +
    # RidgePen) coef
    object$S[[1]][cbind(1:(object$m[1] + 2), 1:(object$m[1] + 2))] <- object$S[[1]][cbind(1:(object$m[1] + 2),
                                                                                          1:(object$m[1] + 2))] + 1
    object$rank <- min(object$bs.dim, object$rank + object$m[1] + 2)
  }
  if (!is.null(object$xt$ridge) && object$xt$ridge) {
    # add lag-lead penalty to first and last <order of B-spline>+1 (=m+2) basis functions penalty
    # = coef' (lambda_1*DiffPen + lambda_2*LagLeadPen) coef
    if (!is.null(object$xt$leadpen) && !is.null(object$xt$lagpen) && object$xt$leadpen && object$xt$lagpen) {
      object$S[[2]] <- matrix(0, object$bs.dim, object$bs.dim)
      # penalize lead
      object$S[[2]][cbind(1:(object$m[1] + 2), 1:(object$m[1] + 2))] <- 1
      # penalize lag
      object$S[[2]][cbind((object$bs.dim - (object$m[1] + 2) + 1): object$bs.dim,
                          (object$bs.dim - (object$m[1] + 2) + 1): object$bs.dim)] <- 1

      object$rank <- c(object$rank, min(2*(object$m[1] + 2), object$bs.dim))
    }
    else if (!is.null(object$xt$leadpen) && object$xt$leadpen) {
      # add ridge penalty (lead penalty) to first <order of B-spline>+1 (=m+2) basis functions penalty
      # = coef' (lambda_1*DiffPen + lambda_2*RidgePen) coef
      object$S[[2]] <- matrix(0, object$bs.dim, object$bs.dim)
      object$S[[2]][cbind(1:(object$m[1] + 2), 1:(object$m[1] + 2))] <- 1
      object$rank <- c(object$rank, object$m[1] + 2)
    }
    else if (!is.null(object$xt$lagpen) && object$xt$lagpen) {
      # add lag penalty to last <order of B-spline>+1 (=m+2) basis functions penalty
      # = coef' (lambda_1*DiffPen + lambda_2*LagPen) coef
      object$S[[2]] <- matrix(0, object$bs.dim, object$bs.dim)
      object$S[[2]][cbind((object$bs.dim - (object$m[1] + 2) + 1): object$bs.dim,
                          (object$bs.dim - (object$m[1] + 2) + 1): object$bs.dim)] <- 1
      object$rank <- c(object$rank, object$m[1] + 2)
    }
  }
  if (!is.null(object$xt$constrain) && object$xt$constrain) {
    # constrain to end in zero (i.e (X%*%coefficients)[1] == 0) -->
    # Constraint matric C = X[1,]
    object$C <- matrix(object$X[1, ], nrow = 1)
    object$C <- structure(object$C, always.apply = TRUE)
  }

  return(object)

}