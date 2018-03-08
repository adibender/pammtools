context("Transformation with TDC")

test_that("split_tdc works correctly", {
	data("pbc", package = "survival")

	event_df <- pbc %>%
	  filter(id <= 5) %>%
	  mutate(event = 1L*(status==2)) %>%
	  select(id, time, event, sex, bili)
	tdc_df <- pbcseq %>%
		filter(id <= 5) %>%
	  select(id, day, bili)

	pbc_ped <- split_tdc(Surv(time, event)~., event_df, tdc_df, te_var="day",
		status_var = "event")


	expect_data_frame(pbc_ped, nrows=93L, ncols=8L)
	expect_is(pbc_ped, "ped")
	expect_subset(c("cut", "id_var", "intvars"), names(attributes(pbc_ped)))

})
