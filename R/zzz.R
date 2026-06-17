.onLoad <- function(
  libname = find.package("pammtools"),
  pkgname = "pammtools"
) {
  if (requireNamespace("mgcv", quietly = TRUE)) {
    base::registerS3method(
      "smooth.construct",
      "fdl.smooth.spec",
      smooth.construct.fdl.smooth.spec,
      envir = asNamespace("mgcv")
    )
  }

  if (getRversion() >= "2.5.1") {
    utils::globalVariables(c(
      ".", "id",
      # NSE column names used in gg_state_occupation()
      "df_long", "prob", "state", "time", "trans_prob_matrix"
    ))
  }
  invisible()
}

.onAttach <- function(
  libname = find.package("pammtools"),
  pkgname = "pammtools"
) {
  invisible()
}
