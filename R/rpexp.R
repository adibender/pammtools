#' Draw random numbers from piece-wise exponential distribution.
#'
#' This is a copy of the same function from \code{rpexp} from package
#' \pkg{msm}.
#' Copied here to reduce dependencies.
#'
#' @inheritParams msm::rpexp
#' @importFrom stats rexp
#'
#' @keywords internal
rpexp <- function (n = 1, rate = 1, t = 0)
{
    if (length(t) != length(rate))
        stop("length of t must be equal to length of rate")
    if (!isTRUE(all.equal(0, t[1])))
        stop("first element of t should be 0")
    if (is.unsorted(t))
        stop("t should be in increasing order")
    if (length(n) > 1)
        n <- length(n)
    if (n == 0)
        return(numeric(0))
    if (length(rate) == 1)
        return(rexp(n, rate))
    ret <- numeric(n)
    left <- 1:n
    for (i in seq_along(rate)) {
        re <- rexp(length(left), rate[i])
        r <- t[i] + re
        success <- if (i == length(rate))
            seq_along(left)
        else which(r < t[i + 1])
        ret[left[success]] <- r[success]
        left <- setdiff(left, left[success])
        if (length(left) == 0)
            break
    }
    ret
}
