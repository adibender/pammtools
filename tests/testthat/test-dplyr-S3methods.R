context("Dplyr methods for specific classes")

data("leuk2", package = "bpcp")
leuk2 <- dplyr::slice(leuk2, 1:2)
ped <- split_data(Surv(time, status)~., data = leuk2, id = "id")

test_that("ped class is preserved", {
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
