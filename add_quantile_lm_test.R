tb <- cars
fit <- lm(cars$dist~cars$speed)

add_quantile.lm(tb,fit, prob = 0.3) %>%
    add_quantile.lm(fit, prob = 0.4) %>%
    head
