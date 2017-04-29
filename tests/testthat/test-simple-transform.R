context("Simple transformation to PED data")

data("leuk2", package="bpcp")
leuk2 <- leuk2 <- leuk2[c(2, 16, 17, 20, 39), ]

test_that("Output as expected without id", {
	expect_is(leuk.ped <- split_data(Surv(time, status)~., data=leuk2, 
		cut=c(0:5, 10, 40)), "ped")
	expect_equal(dim(leuk.ped), c(27, 7))
	expect_is(leuk.ped, "data.frame")
	expect_true(all(c("ped_status", "tstart", "tend", "interval", "offset") %in% 
		names(leuk.ped)))
	expect_is(attr(leuk.ped, "cut"), "numeric")
	expect_is(attr(leuk.ped, "intvars"), "character")
})



test_that("Output as expected with id", {
	expect_is(leuk.ped <- split_data(Surv(time, status)~., data=leuk2, 
		cut=c(0:5, 10, 40), id="id"), "ped")
	expect_equal(dim(leuk.ped),c(27, 8))
})

test_that("ID kept when id variable excluded in formula", {
	expect_is(leuk.ped <- split_data(Surv(time, status)~treatment, data=leuk2, 
		cut=c(0:5, 10, 40), id="id"), "ped")
	expect_equal(ncol(leuk.ped), 7)
})


test_that("Error on wrong input", {
	expect_error(split_data())
	expect_error(split_data(x~y, data=leuk2, cut=c(0:5, 10, 40)))
	expect_error(split_data(Surv(time2, status)~., data=leuk2, cut=c(0:5, 10, 40)))
	expect_error(split_data(Surv(ped_time, status)~., 
		data=rename(leuk2, ped_time=time))) # already in data set ped_time
})


# str(leuk.ped)
# leuk.ped
# cph <- coxph(Surv(time, status)~treatment, data=leuk2)
# pem <- glm(status ~ -1 + interval + treatment, data=leuk.ped, offset=offset, 
# 	family=poisson())
# summary(pem)
# coef(cph)
# coef(pem)["treatmentplacebo"]
# # summary(cph.ped <- coxph(Surv(time, status)~treatment, data=leuk.ped))

# cut=unique(leuk2$time)
# leuk.ped <- split_data(Surv(time, status)~., data=leuk2, cut=cut)
# pem <- glm(status ~ interval + treatment, data=leuk.ped, offset=offset, family=poisson())
# coef(cph)
# summary(pem)


# # ## compare baselines 
# cph.base <- basehaz(cph, center=FALSE)
# plot(cph.base$time, cph.base$hazard, type="s")
# b <- coef(pem)
# h <- exp(b[1] + c(0, b[2:24]))
# H <- cumsum(h * diff(c(0, cph.base$time)))
# S <- exp(-H)
# lines(cph.base$time, H, col=2, type="s")

# ndf <- data.frame(interval=unique(leuk.ped$interval), treatment="placebo")
# ndf2 <- data.frame(interval=unique(leuk.ped$interval), treatment="6-MP")

# bp <- predict(pem, newdata=ndf, type="response")
# bp2 <- predict(pem, newdata=ndf2, type="response")
