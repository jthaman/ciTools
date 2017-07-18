rm(list =ls())
library(lme4)
library(ciTools)
library(tidyverse)

set.seed(20170718)
source("simulations/sim functions.R")

system.time(results <- sim_study_ci(n = 250, fixed_ef = c(1,1), sigma_y = 5,
                                    sigma_g = 1, group_no = 5,
                                    sims = 10, bootSims = 200,
                                    includeRanef = TRUE, alpha = 0.5));results

