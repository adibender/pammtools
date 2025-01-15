#' pammtools: Piece-wise exponential Additive Mixed Modeling tools.
#'
#' \code{pammtools} provides functions and utilities that facilitate fitting
#' Piece-wise Exponential Additive Mixed Models (PAMMs), including data
#' transformation and other convenience functions for pre- and post-processing
#' as well as plotting.
#'
#' The best way to get an overview of the functionality provided and how to
#' fit PAMMs is to view the vignettes
#' available at \url{https://adibender.github.io/pammtools/articles/}.
#' A summary of the vignettes' content is given below:
#'
#' \itemize{
#' \item \href{https://adibender.github.io/pammtools/articles/basics.html}{basics}:
#' Introduction to PAMMs and basic modeling.
#' \item \href{https://adibender.github.io/pammtools/articles/baseline.html}{baseline}:
#' Shows how to estimate and visualize baseline model (without covariates) and
#' comparison to respective Cox-PH model.
#' \item \href{https://adibender.github.io/pammtools/articles/convenience.html}{convenience}:
#' Convenience functions for post-processing and plotting PAMMs.
#' \item \href{https://adibender.github.io/pammtools/articles/data-transformation.html}{data-transformation}:
#' Transforming data into a format suitable to fit PAMMs.
#' \item \href{https://adibender.github.io/pammtools/articles/frailty.html}{frailty}:
#' Specifying "frailty" terms, i.e., random effects for PAMMs.
#' \item \href{https://adibender.github.io/pammtools/articles/splines.html}{splines}:
#' Specifying spline smooth terms for PAMMs.
#' \item \href{https://adibender.github.io/pammtools/articles/strata.html}{strata}:
#' Specifying stratified models in which each level of a grouping variable has a
#' different baseline hazard.
#' \item \href{https://adibender.github.io/pammtools/articles/tdcovar.html}{tdcovar}:
#' Dealing with time-dependent covariates.
#' \item \href{https://adibender.github.io/pammtools/articles/tveffects.html}{tveffects}:
#' Specifying time-varying effects.
#' \item \href{https://adibender.github.io/pammtools/articles/left-truncation.html}{left-truncation}:
#' Estimation for left-truncated data.
#'\item \href{https://adibender.github.io/pammtools/articles/competing-risks.html}{competing-risks}:
#' Competing risks analysis.
#' }
#'
#' @name pammtools
#' @keywords internal
"_PACKAGE"
#' @references
#' Bender, Andreas, Andreas Groll, and Fabian Scheipl. 2018.
#' “A Generalized Additive Model Approach to Time-to-Event Analysis”
#' Statistical Modelling, February. https://doi.org/10.1177/1471082X17748083.
#'
#' Bender, Andreas, Fabian Scheipl, Wolfgang Hartl, Andrew G. Day, and Helmut Küchenhoff. 2019.
#' “Penalized Estimation of Complex, Non-Linear Exposure-Lag-Response Associations.”
#' Biostatistics 20 (2): 315–31. https://doi.org/10.1093/biostatistics/kxy003.
#'
#' Bender, Andreas, and Fabian Scheipl. 2018.
#' “pammtools: Piece-Wise Exponential Additive Mixed Modeling Tools.”
#' ArXiv:1806.01042 [Stat], June. https://arxiv.org/abs/1806.01042.
NULL
NULL
