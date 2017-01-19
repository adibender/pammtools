.onLoad <- function(libname=find.package("pam"), pkgname="pam") {
	
	if(getRversion() >= "2.5.1") {
		utils::globalVariables(
			c("status", "time", "tstart"))
	}

	invisible()

}