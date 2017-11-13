.onLoad <- function(libname=find.package("pammtools"), pkgname="pammtools") {

	if(getRversion() >= "2.5.1") {
		utils::globalVariables(
			c("tstart", "tend", "interval", "intlen", "intmid",
				".", "eff", "high", "low", "type", "se", "tmp.fit", "tmp.se", "hazard",
				"lower", "upper", "term", "x", "fit", "se.fit",
				"ped_censored", "idx_train", "idx_test"))
	}

	invisible()

}
