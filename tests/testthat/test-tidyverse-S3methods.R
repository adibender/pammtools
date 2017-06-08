context("Tidyverse methods for specific classes")

test_that("ped class is preserved after dplyr operations", {
	data("veteran", package = "survival")
	veteran <- dplyr::slice(veteran, 1:2)
	ped <- split_data(Surv(time, status)~ trt + age, cut=c(0, 100, 400), 
		data = veteran, id = "id")

	expect_is(filter(ped, id == 1), "ped")
	expect_is(slice(ped, 1), "ped")
	expect_is(arrange(ped, desc(id)), "ped")
	expect_is(select(ped, id), "ped")
	expect_is(rename(ped, ID=id), "ped")
	expect_is(mutate(ped, id = id + 1), "ped")
	expect_is(transmute(ped, id = id + 1), "ped")
	expect_is(sample_n(ped, 1), "ped")
	expect_is(sample_frac(ped, 0.5), "ped")
})



test_that("ped class is preserved after tidyr operations", {

	data("pbc", package = "survival")
	tdc_df <- pbcseq %>% 
		filter(id %in% 1:5) %>% 
	  select(id, day, chol)
	class(tdc_df) <- c("ped", class(tdc_df))

	## standard verstion
	temp <- tidyr::fill(tdc_df, chol)
	expect_is(temp, "ped")
	expect_identical(sum(is.na(temp)), 0L)

	## non standard evaluation
	temp <- tidyr::fill_(tdc_df, c("chol"))
	expect_is(temp, "ped")
	expect_identical(sum(is.na(temp)), 0L)

})

