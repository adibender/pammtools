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
    utils::globalVariables(c(".", "id"))
  }
  invisible()
}

.onAttach <- function(
  libname = find.package("pammtools"),
  pkgname = "pammtools"
) {
  invisible()
}
