.onLoad <- function(libname=find.package("pammtools"), pkgname="pammtools") {

  if (getRversion() >= "2.5.1") {

    utils::globalVariables(".")

  }
  invisible()

}

.onAttach <- function(libname=find.package("pammtools"), pkgname="pammtools") {

  invisible()
}
