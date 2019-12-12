hweibull <- function(t, alpha=1.5, tau=2) {
  (alpha/tau) * (t/tau)^(alpha -1)
}

sdata <- data.frame(time=rweibull(500, 1.5, 2))
sdata$event = 1

kappa5 <- seq(0, 4, by=0.8)
base_ped <- pammtools::split_data(Surv(time, event)~., data=sdata, cut=kappa5,
 id="id")
usethis::use_data(base_ped)
