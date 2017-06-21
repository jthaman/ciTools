library(tidyverse)

tb <- cars
fit <- lm(cars$dist~cars$speed)

tb <- add_quantile.lm(tb,fit, prob = 0.025) %>%
    add_quantile.lm(fit, prob = 0.975)

tb <- add_quantile.lm(tb,fit, prob = 0.025)

tb <- add_pi.lm(tb, fit)

tb %>% head ## the same

ggplot(data = tb, aes(x = speed, y = dist)) +
    geom_point() +
    geom_line(data = tb, aes(x = speed, y = pred)) +
    geom_line(data = tb, aes(x = speed, y = quantile0.025), color = "blue") +
    geom_line(data = tb, aes(y = quantile0.4), color = "red") +
    geom_line(data = tb, aes(y = quantile0.975), color = "green") 


## ggplot example
ggplot(tb, aes(x = speed)) +
   geom_ribbon(aes(ymin = LPB, ymax = UPB),
               fill = "blue", alpha = 0.2) +
   geom_point(aes(y = dist)) +
   geom_line(aes(y = pred), colour = "blue", size = 1) 
