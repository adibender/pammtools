staph <- readxl::read_excel("staph.xlsx") %>%
  mutate(id = as.factor(id))

usethis::use_data(staph, overwrite = TRUE)
