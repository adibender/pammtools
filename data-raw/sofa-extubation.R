library(dplyr)
library(tidyr)
data(VAP_data, package = "TBFmultinomial")
# create unique IDs
VAP_data <- VAP_data %>%
  mutate(tmp_id =  paste(ID, day, sep = ".")) %>%
  group_by(ID) %>%
  mutate(csdup = cumsum(duplicated(tmp_id))) %>%
  ungroup() %>%
  mutate(ID = ifelse(csdup > 0, ID + 1000, ID))

# assume constant SOFA score between updates
VAP_complete <- VAP_data %>%
  group_by(ID) %>%
  mutate(time = max(day)) %>%
  ungroup() %>%
  complete(ID, day = full_seq(day, 1)) %>%
  fill(gender, type, SAPSadmission, SOFA, outcome, time, .direction = "down") %>%
  filter(day <= time)

extub_event <- VAP_complete %>%
  select(ID, gender, type, SAPSadmission, outcome, time) %>%
  group_by(ID) %>%
  slice(n()) %>%
  mutate(extubation = 1 * (outcome == "extubated")) %>%
  ungroup() %>%
  select(-outcome)

extub_tdc <- VAP_complete %>%
  select(ID, day, SOFA) %>%
  mutate(day = day - 1) %>% # assume that SOFA available at the beginning of the day
  filter(day <= 49)

use_data(extub_event, extub_tdc, overwrite = TRUE)
