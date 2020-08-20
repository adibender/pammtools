data("sir.adm", package = "mvna")

sir_adm <- sir.adm[c(1,
  2, 3, 26, 40, 43, 50), ]


sir_adm2 <- sir.adm[1:150, ]

usethis::use_data(sir_adm, sir_adm2, internal = TRUE, overwrite = TRUE)
