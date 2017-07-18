library(lme4)
library(ciTools)
library(tidyverse)

## simulation study comparison of LMM CI techniques in ciTools

sim_study_ci <- function(n, fixed_ef, sigma_y, sigma_g, group_no, sims,
                         bootSims, includeRanef, alpha){

    success_para <- success_pi <- success_boot <- rep(0, sims)
    ## model params
    b0 <- fixed_ef[1]
    b1 <- fixed_ef[2]
    group <- rep (1:group_no, each = n/group_no) 
    a <- rnorm (group_no, 0, sigma_g) 
    x <- seq(1,2, length.out = n)
    EYx <- b0 + b1*x 
    EYxg <- b0 + b1*x + a[group]

    if (includeRanef)
        truth <- EYxg
    else
        truth <- EYx
    
    for (i in 1:sims){
        y <- rnorm (n, fixed_ef[1] + a[group] + fixed_ef[2]*x, sigma_y) 
        tb <<- as_data_frame(cbind (x, group, y)) ## I hate this
        
        fit <- lmer(y ~ x + (1 | group))

        tb <- add_ci(tb, fit, includeRanef = includeRanef, alpha = alpha, names = c("lwr_para", "upr_para")) %>%
            add_ci(fit, type = "sim", includeRanef = includeRanef, alpha = alpha, nSims = 1000, names = c("lwr_pi", "upr_pi")) %>%
            add_ci(fit, type = "boot", includeRanef = includeRanef, alpha = alpha, nSims = bootSims, names = c("lwr_boot", "upr_boot"))

        if (all( (tb$lwr_para < truth) && (tb$upr_para > truth) ))
            success_para[i] <- 1 
        if (all( (tb$lwr_pi < truth) && (tb$upr_pi > truth) ))
            success_pi[i] <- 1 
        if (all( (tb$lwr_boot < truth) && (tb$upr_boot > truth) ))
            success_boot[i] <- 1 
        rm(tb)

    }
    return(list(para_coverage_prob = mean(success_para),
                pi_coverage_prob = mean(success_pi),
                boot_coverage_prob = mean(success_boot)))
}



system.time(results <- sim_study_ci(n = 500, fixed_ef = c(1,1), sigma_y = 4,
                         sigma_g = 1, group_no = 5,
                         sims = 100, bootSims = 200,
                         includeRanef = TRUE, alpha = 0.5)); results

## n = 500 and sims = 25 : 115 secs
## n = 500 and sims = 100 : 8 mins

sim_study_pi <- function(n, fixed_ef, sigma_y, sigma_g, group_no, sims,
                         bootSims, includeRanef, alpha){

    success_para <- success_pi <- success_boot <- rep(0, sims)
    ## model params
    b0 <- fixed_ef[1]
    b1 <- fixed_ef[2]
    group <- rep (1:group_no, each = n/group_no) 
    a <- rnorm (group_no, 0, sigma_g) 
    x <- seq(1,2, length.out = n)
    EYx <- b0 + b1*x 
    EYxg <- b0 + b1*x + a[group]

    if (includeRanef)
        truth <- EYxg
    else
        truth <- EYx
    
    for (i in 1:sims){
        y <- rnorm (n, fixed_ef[1] + a[group] + fixed_ef[2]*x, sigma_y) 
        tb <- as_data_frame(cbind (x, group, y)) ## I hate this
        
        fit <- lmer(y ~ x + (1 | group))

        tb <- add_pi(tb, fit, includeRanef = includeRanef, alpha = alpha, names = c("lwr_para", "upr_para")) %>%
            add_pi(fit, type = "sim", includeRanef = includeRanef, alpha = alpha, nSims = 1000, names = c("lwr_pi", "upr_pi")) %>%
            add_pi(fit, type = "boot", includeRanef = includeRanef, alpha = alpha, nSims = bootSims, names = c("lwr_boot", "upr_boot"))

        if (all( (tb$lwr_para < truth) && (tb$upr_para > truth) ))
            success_para[i] <- 1 
        if (all( (tb$lwr_pi < truth) && (tb$upr_pi > truth) ))
            success_pi[i] <- 1 
        if (all( (tb$lwr_boot < truth) && (tb$upr_boot > truth) ))
            success_boot[i] <- 1 
        rm(tb)

    }
    return(list(para_coverage_prob = mean(success_para),
                pi_coverage_prob = mean(success_pi),
                boot_coverage_prob = mean(success_boot)))
}








