context("Transformation with TDC")

test_that("split_tdc works correctly", {
	data("pbc", package = "survival")

	event_df <- pbc %>% 
	  filter(id <= 5) %>% 
	  mutate(event = 1*(status==2)) %>% 
	  select(id, time, event, sex, bili)

## we rename the "day" variable here for use of `split_tdc` function later
	tdc_df <- pbcseq %>% 
		filter(id %in% 1:5) %>% 
	  select(id, day, bili) %>% 
	  rename(time = day)
	pbc_ped <- split_tdc(Surv(time, event)~., event_df, tdc_df, 
		"id", "time", "status")

	expect_data_frame(pbc_ped, nrows=85L, ncols=8L)
	expect_is(pbc_ped, "ped")
	expect_subset(c("cut", "id_var", "intvars"), names(attributes(pbc_ped)))

})