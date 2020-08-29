library(pammtools)
data(tumor)
ped_tumor <- tumor %>% as_ped(Surv(days, status)~., max_time =3034)
pam_tumor_tve <- mgcv::bam(
  formula = ped_status ~ ti(tend) +
    complications + ti(tend, by = as.ordered(complications)) +
    metastases + ti(tend, by = as.ordered(metastases)) +
    sex + ti(tend, by = as.ordered(sex)) +
    transfusion + ti(tend, by = as.ordered(transfusion)) +
    resection + ti(tend, by = as.ordered(resection)) +
    s(tend, by = charlson_score) +
    s(tend, by = age),
  data = ped_tumor, family = poisson(), offset = offset,
  method = "fREML", discrete = TRUE)

cumu_coefs <- get_cumu_coef(pam_tumor_tve, ped_tumor, terms = c("metastases", "charlson_score", "complications", "sex"))

library(timereg)
aal <- aalen(Surv(days, status)~., data = tumor)
cumu_coef_aalen <- get_cumu_coef(aal, terms = c("metastases", "charlson_score", "complications", "sex"))
cumu_coef_aalen$variable[cumu_coef_aalen$variable == "complicationsyes"] <- "complications (yes)"
cumu_coef_aalen$variable[cumu_coef_aalen$variable == "metastasesyes"] <- "metastases (yes)"
cumu_coef_aalen$variable[cumu_coef_aalen$variable == "sexfemale"] <- "sex (female)"

c1 <- cumu_coefs %>% filter(variable == "complications (yes)")
c2 <- cumu_coef_aalen %>% filter(variable == "complications (yes)")

library(ggplot2)
theme_set(theme_minimal())
p <- ggplot(c1, aes(x = time, y = cumu_hazard)) +
  geom_hazard(data = c2, aes(y = cumu_hazard), col = "grey30", alpha = .7) +
  geom_ribbon(data = c2, aes(ymin = cumu_upper, ymax = cumu_lower), alpha = .3, fill = "grey70") +
  geom_hazard(col = "steelblue", alpha = .7) +
  geom_ribbon(aes(ymin = cumu_lower, ymax = cumu_upper), alpha = .3, fill = "steelblue") +
  ylab("") + xlab("") + theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_blank()
  ) +
  facet_wrap(~variable, nrow = 2)

library(hexSticker)
sticker(p, package = "pammtools", p_size = 20, s_width = 2, s_height = 2,
  p_color = "black", h_fill = "whitesmoke", h_color = "black",
  s_x = .9, s_y = 1)
