#' Calculate difference in cumulative hazards and respective standard errors
#'
#' Uses Delta method to calculate standard errors for cumulative hazard
#' differences based on coefficients and standard errors of log hazard
#' (code adapted from fabian-s used in consulting project).
#'
#' @param d1 A data set used as \code{newdata} in \code{predict.gam}
#' @param d2 See \code{d1}
#' @param model A model object for which a predict method is implemented which
#' returns the design matrix (e.g., \code{mgcv::gam}).
#' @importFrom mgcv predict.gam
#' @keywords internal
compute_cumu_diff <-  function(d1, d2, model) {
  X1 <- predict.gam(model, newdata = d1, type = "lpmatrix")
  X2 <- predict.gam(model, newdata = d2, type = "lpmatrix")
  haz1 <- exp(drop(X1 %*% model$coefficients))
  haz2 <- exp(drop(X2 %*% model$coefficients))
  cumu_diff <- cumsum(haz2*d2$intlen) - cumsum(haz1*d1$intlen)
  # delta rule:
  # transformation:
  # f(t, coef) = S(t|1)/S(t|2) =
  #         = exp(- int^t_0 lambda(s|1)ds + int^t_0 lambda(s|2)ds)
  #         = exp((1,...,1,0,...,0) %*% (exp(X2%*%coef) - exp(X1%*%coef))
  ## (pre-multiplication with 0-1-vec is (cum)sum)
  # d/d coef f(t, coef) = f(t, coef) * (X2'((1,...,1,0,...,0)*exp(X2%*%coef)) -
  #                                     X1'((1,...,1,0,...,0)*exp(X1%*%coef)))
  jacobi <- matrix(NA, ncol=length(model$coefficients), nrow = length(cumu_diff))
  lo.tri <- lower.tri(diag(nrow(X1)), diag=TRUE)
  for(day in 1:nrow(jacobi)) {
    jacobi[day,] <- t(X2)%*%(lo.tri[day,]*haz2*d2$intlen) - t(X1)%*%(lo.tri[day,]*haz1*d1$intlen) *
      cumu_diff[day]
    }
  se_cdiff <- sqrt(diag(jacobi %*% model$Vp %*% t(jacobi)))
  list(cumu_diff = cumu_diff, se_cumu = se_cdiff)

  }
