tb <- cars
fit <- lm(dist ~ speed, data = cars)
add_ci(tb, fit)
