dat <- lme4::sleepstudy
fit1 <- lm(Reaction ~ Days, data = dat)
?add_ci
add_ci(dat, fit1) %>% add_ci(fit1)
?add_pi
add_pi(dat, fit1) %>% add_pi(fit1)
?add_probs
?add_quantile
add_quantile(dat, fit1, p = 0.6) %>%
    add_quantile(fit1, p = 0.7)

add_ci(dat, fit1) %>%
  add_ci(fit1, alpha = 0.1) %>%
  add_ci(fit1, alpha = 0.5)

# Linear Regression with Missing Data
dat_missing <- sleepstudy
dat_missing[3,1] <- NA
dat_missing[4,2] <- NA
add_ci(dat_missing, fit1)

# Linear Mixed Models
library(lme4)
dat <- sleepstudy
fit2 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
add_ci(dat, fit2)
add_ci(dat, fit2, type = "sim")
add_pi(dat, fit2)
add_pi(dat, fit2, type = "sim") %>%
add_pi(fit2, type = "parametric")
add_probs(dat, fit2, q = 400, comparison = ">")
add_probs(dat, fit2, q = 400, comparison = ">", type = "sim", nSims = 1000) %>%
add_probs(fit2, q = 400, comparison = ">", type = "boot", nSims = 1000)
add_quantile(dat, fit2, p = 0.7)
add_quantile(dat, fit2, p = 0.7, type = "boot")
add_quantile(dat, fit2, p = 0.7, type = "sim")

# Without Random Effects
add_ci(dat, fit2, includeRanef = FALSE)
add_ci(dat, fit2, type = "sim", includeRanef = FALSE)
add_pi(dat, fit2, includeRanef = FALSE)
add_pi(dat, fit2, type = "sim", includeRanef = FALSE)

#logistic Regression
dat2 <- cars
fit3 <- glm(I(dist > 20) ~ speed, data = cars, family = "binomial")
add_ci(dat2, fit3)
add_pi(dat2, fit3)
add_probs(dat2, fit3, q = 0, comparison = "=")
add_quantile(dat2, fit3, p = 0.7)

#logistic regression warning
fit3 <- glm(I(dist > 2) ~ speed, data = cars, family = "binomial")
add_ci(dat2, fit3)

#Poisson Regression
fit4 <- glm(dist ~ speed, data = cars, family = "poisson")
add_ci(dat2, fit4)
add_pi(dat2, fit4)
add_probs(dat2, fit4, q = 50, comparison = "<")
add_quantile(dat2, fit4, p = 0.7)

#Log-Linear Regession
fit5 <- lm(log(exp(dist)) ~ speed, data = cars)
add_ci(dat2, fit5, log_response = TRUE)
add_pi(dat2, fit5, log_response = TRUE)
add_probs(dat2, fit5, q = 50, log_response = TRUE)
add_quantile(dat2, fit5, p = 0.6, log_response = TRUE)
