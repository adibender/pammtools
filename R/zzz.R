
.onLoad <- function(libname=find.package("pammtools"), pkgname="pammtools") {

  if (getRversion() >= "2.5.1") {

    utils::globalVariables(".")

  }

  invisible()

}

.onAttach <- function(libname=find.package("pammtools"), pkgname="pammtools") {

  packageStartupMessage(
    "The add_term function has changed in pammtools 0.1.12.
    The argument `relative` is deprecated and replaced by `reference`.
    See help page for details.")

  invisible()
}
