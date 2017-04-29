.onLoad <- function(libname=find.package("pam"), pkgname="pam") {
	
	if(getRversion() >= "2.5.1") {
		utils::globalVariables(
			c("ped_status", "ped_time", "tstart", "tend", "interval", "intlen", "intmid", 
				".", "eff", "se", "tmp.fit", "tmp.se", "hazard", "lower", "upper", "term",
				"x", "ci.lower", "ci.upper", "fit", "se.fit"))
	}

	invisible()

}