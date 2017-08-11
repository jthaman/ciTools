## ------------------------------------------------------------------------
## render("u:/org/ciTools-demo.R")
## setup 
rm(list = ls())
set.seed(20170808)
knitr::opts_chunk$set(fig.width=6, fig.height=4) 
library(rmarkdown)
library(ciTools)
library(tidyverse)
library(knitr)

## ------------------------------------------------------------------------
my_data <- cars
glimpse(my_data)

## ------------------------------------------------------------------------
model <- lm(dist ~ speed, data = cars)

## ------------------------------------------------------------------------
my_data_with_ci <- add_ci(my_data, model, names = c("lcb", "ucb"))
kable(head(my_data_with_ci, n =10), row.names = TRUE)

## ------------------------------------------------------------------------
my_data_with_ci %>%
    ggplot(aes(x = speed, y = dist)) +
    geom_point(size = 2) +
    geom_line(aes(y = pred), size = 2, color = "maroon") +
    geom_ribbon(aes(ymin = lcb, ymax = ucb), fill =
    "royalblue1", alpha = 0.3) + 
    ggtitle("Stopping Distance vs. Car Speed: 95% Confidence Interval") +
    xlab("Car Speed (mph)") +
    ylab("Stopping Distance (ft)")

## ------------------------------------------------------------------------
my_data_with_pi <- add_pi(my_data, model, names = c("lpb", "upb"))

## ------------------------------------------------------------------------
kable(head(my_data_with_pi, n = 10), row.names = TRUE)

## ------------------------------------------------------------------------
my_data %>%
  add_ci(model, names = c("lcb", "ucb")) %>%
  add_pi(model, names = c("lpb", "upb")) %>%
    ggplot(aes(x = speed, y = dist)) +
    geom_point(size = 2) +
    geom_line(aes(y = pred), size = 2, color = "maroon") +
    geom_ribbon(aes(ymin = lpb, ymax = upb), fill = "orange2",
                alpha = 0.3) +
    geom_ribbon(aes(ymin = lcb, ymax = ucb), fill =
    "royalblue1", alpha = 0.3) + 
    ggtitle("Stopping Distance vs. Car Speed: 95% CI and 95% PI") +
    xlab("Car Speed (mph)") +
    ylab("Stopping Distance (ft)")

## ------------------------------------------------------------------------
my_data %>%
  add_probs(model, q = 70) %>%
    ggplot(aes(x = speed, y = prob_less_than70)) +
    geom_line(aes(y = prob_less_than70), size = 2, color = "maroon") +
    scale_y_continuous(limits = c(0,1)) +
    ggtitle("Probability Stopping Distance is Less Than 70") +
    xlab("Car Speed (mph)") +
    ylab("Pr(Dist < 70)")

## ------------------------------------------------------------------------
my_data %>%
  add_pi(model, names = c("lpb", "upb")) %>%
  add_quantile(model, p = 0.9) %>%
    ggplot(aes(x = speed, y = dist)) +
    geom_point(size = 2) +
    geom_line(aes(y = pred), size = 2, color = "maroon") +
    geom_line(aes(y = quantile0.9), size = 2, color = "forestgreen") + 
    geom_ribbon(aes(ymin = lpb, ymax = upb), fill = "orange2",
                alpha = 0.3) +
    ggtitle("Stopping Distance vs. Car Speed: 95% PI with 0.9-Quantile") +
    xlab("Car Speed (mph)") +
    ylab("Stopping Distance (ft)")

## ------------------------------------------------------------------------
data_with_results <- my_data %>%
  add_ci(model) %>%
  add_pi(model) %>%
  add_probs(model, q= 70) %>%
  add_quantile(model, p = 0.9) 

kable(head(data_with_results))

