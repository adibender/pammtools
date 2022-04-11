#' New basis for penalized lag selection

smooth.construct.fdl.smooth.spec <- function(object, data, knots) {

  # modify object so that it's fitted as a p-spline signal regression term:
  object$bs <- "ps"
  object <- mgcv:::smooth.construct.ps.smooth.spec(object, data, knots)

  if (!is.null(object$xt$fullrankpen) && object$xt$fullrankpen) {
    # add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions with
    # same variance as difference penalty : penalty = lambda * coef' (DiffPen +
    # RidgePen) coef
    object$S[[1]][cbind(1:(object$m[1] + 2), 1:(object$m[1] + 2))] <- object$S[[1]][cbind(1:(object$m[1] +
      2), 1:(object$m[1] + 2))] + 1
    object$rank <- min(object$bs.dim, object$rank + object$m[1] + 2)
  }
  if (!is.null(object$xt$ridge) && object$xt$ridge) {
    # add ridge penalty to first <order of B-spline>+1 (=m+2) basis functions penalty
    # = coef' (lambda_1*DiffPen + lambda_2*RidgePen) coef
    object$S[[2]] <- matrix(0, object$bs.dim, object$bs.dim)
    object$S[[2]][cbind(1:(object$m[1] + 2), 1:(object$m[1] + 2))] <- 1
    object$rank <- c(object$rank, object$m[1] + 2)
  }
  if (!is.null(object$xt$constrain) && object$xt$constrain) {
    # constrain to end in zero (i.e (X%*%coefficients)[1] == 0) -->
    # Constraint matric C = X[1,]
    object$C <- matrix(object$X[1, ], nrow = 1)
    object$C <- structure(object$C, always.apply = TRUE)
  }

  return(object)

}
