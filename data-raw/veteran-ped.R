hweibull <- function(t, alpha=1.5, tau=2) {
  (alpha/tau) * (t/tau)^(alpha -1)
}

tdf <- data.frame(time = seq(0, 4, by=0.05))
tdf$hazard 
	mutate(hweib = hweibull(time))

set.seed(16042017)
sdata <- data.frame(time=rweibull(500, 1.5, 2))
sdata$event = 1

kappa5 <- seq(0, 4, by=0.8)
base_ped <- split_data(Surv(time, event)~., data=sdata, cut=kappa5, id="id")
use_data(base_ped)