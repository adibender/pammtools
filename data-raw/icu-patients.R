library(dplyr)

event_df <- readRDS(url("https://github.com/adibender/elra-biostats/raw/master/data/patient.Rds"))
event_df <- event_df%>%
 select(Year:CombinedID, Survdays, PatientDied, survhosp, Gender, Age,
  AdmCatID, ApacheIIScore, BMI, DiagID2) %>%
 mutate(
  Survdays = round(Survdays, 1),
  survhosp = round(survhosp, 1))

tdc_df <- readRDS(url("https://github.com/adibender/elra-biostats/raw/master/data/mergedAndCleanedData.Rds"))

tdc_df <- tdc_df %>%
  select(CombinedID, Study_Day, caloriesPercentage, proteinGproKG)

RNGkind(sample.kind = "Rounding")
set.seed(5032018)
samp_id <- sample(unique(tdc_df$CombinedID), 2000)

patient <- filter(event_df, CombinedID %in% samp_id) %>%
  arrange(CombinedID)

daily <- filter(tdc_df, CombinedID %in% samp_id) %>%
  arrange(CombinedID)

usethis::use_data(patient, daily, overwrite = TRUE)
