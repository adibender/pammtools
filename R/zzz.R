.onLoad <- function(libname=find.package("pam"), pkgname="pam") {
	
	if(getRversion() >= "2.5.1") {
		utils::globalVariables(
			c("ped_status", "ped_time", "tstart", "tend", "interval", "intlen", "intmid", 
				".", "eff", "high", "low", "type", "se", "tmp.fit", "tmp.se", "hazard", 
				"lower", "upper", "term", "x", "ci.lower", "ci.upper", "fit", "se.fit", 
				"ped_censored", "cumhazard", "cumlower", "cumupper", 
				"idx_train", "idx_test"))
	}

	invisible()

}