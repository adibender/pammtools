# competing risks data
data("sir.adm", package = "mvna")

sir_adm <- sir.adm[c(1, 2, 3, 26, 40, 43, 50), ]


sir_adm2 <- sir.adm[1:150, ]

# usethis::use_data(sir_adm, sir_adm2, internal = TRUE, overwrite = TRUE)

# recurrent events data
data("cgd", package = "frailtyHL")
cgd <- cgd[cgd$enum <= 2, ]
cgd2 <- cgd[cgd$id %in% c(1:2), ]
# usethis::use_data(cgd, cgd2, internal = TRUE, overwrite = TRUE)

# left truncated data
data("mort", package = "eha")
mort <- mort %>% group_by(id) %>% slice(1) %>% rename(tstart = enter) %>%
  mutate(
    tstart = round(tstart, 3),
    exit = round(exit, 3))

mort2 <- mort %>% filter(id %in% c(1:2))

usethis::use_data(sir_adm, sir_adm2, cgd, cgd2, mort, mort2, internal = TRUE, overwrite = TRUE)
