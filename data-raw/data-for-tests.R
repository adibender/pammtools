data("sir.adm", package = "mvna")
sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]

usethis::use_data(sir_adm, internal = TRUE, overwrite = TRUE)
