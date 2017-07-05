library(tidyverse)
library(lme4)
#library(devtools)
#library(roxygen2)
library(addr)
dat <- sleepstudy

# Linear Regression
fit1 <- lm(Reaction ~ Days, data = sleepstudy)
add_ci(dat, fit1)
add_pi(dat, fit1)
add_probs(dat, fit1, quant = 200)
add_quantile(dat, fit1, p = 0.6)

add_ci(dat, fit1) %>%
  add_ci(fit1, alpha = 0.1) %>%
  add_ci(fit1, alpha = 0.5)

# Linear Regression with Missing Data
dat_missing <- sleepstudy
dat_missing[3,1] <- NA
dat_missing[4,2] <- NA
add_ci(dat_missing, fit1)

# Linear Mixed Models
fit2 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
add_ci(dat, fit2)
add_ci(dat, fit2, ciType = "sim")
add_pi(dat, fit2)
add_pi(dat, fit2, piType = "sim")
add_probs(dat, fit2, quant = 400, comparison = ">")
add_quantile(dat, fit2, prob = 0.7)

# Without Random Effects
add_ci(dat, fit2, includeRanef = FALSE)
add_ci(dat, fit2, ciType = "sim", includeRanef = FALSE)
add_pi(dat, fit2, includeRanef = FALSE)
add_pi(dat, fit2, piType = "sim", includeRanef = FALSE)

#logistic Regression
dat2 <- cars
fit3 <- glm(I(dist > 20) ~ speed, data = cars, family = "binomial")
add_ci(dat2, fit3)
add_pi(dat2, fit3)
add_probs(dat2, fit3, quant = 0, comparison = "=")
add_quantile(dat2, fit3, prob = 0.7)

#logistic regression warning
fit3 <- glm(I(dist > 2) ~ speed, data = cars, family = "binomial")
add_ci(dat2, fit3)

#Poisson Regression
fit4 <- glm(dist ~ speed, data = cars, family = "poisson")
add_ci(dat2, fit4)
add_pi(dat2, fit4)
add_probs(dat2, fit4, quant = 50, comparison = "<")
add_quantile(dat2, fit4, prob = 0.7)

#Log-Linear Regession
fit5 <- lm(log(exp(dist)) ~ speed, data = cars)
add_ci(dat2, fit5, log_response = TRUE)
add_pi(dat2, fit5, log_response = TRUE)
add_probs(dat2, fit5, quant = 50, log_response = TRUE)
add_quantile(dat2, fit5, p = 0.6, log_response = TRUE)
