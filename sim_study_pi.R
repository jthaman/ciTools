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
        y <- rnorm (n, b0 + a[group] + b1*x, sigma_y) 
        tb <<- as_data_frame(cbind (x, group, y)) 
        
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

system.time(results <- sim_study_ci(n = 500, fixed_ef = c(1,1), sigma_y = 1,
                         sigma_g = 1, group_no = 50, sims = 50, bootSims = 200,
                         includeRanef = TRUE, alpha = 0.5)); results

## n = 500, sims = 25, bootsims = 200 : 115 secs
## n = 500 and sims = 100, bootsims = 200 : 8 mins

sim_study_pi <- function(n, fixed_ef, sigma_y, sigma_g, group_no, sims,
                         bootSims, includeRanef, alpha){

    store_para <- store_pi <- store_sim <- 0

    ## model params
    b0 <- fixed_ef[1]
    b1 <- fixed_ef[2]
    group <- rep (1:group_no, each = n/group_no) 
    a <- rnorm (group_no, 0, sigma_g) 
    x <- seq(1,2, length.out = n)
    EYx <- b0 + b1*x 
    EYxg <- b0 + b1*x + a[group]
    y <- rnorm (n, b0 + a[group] + b1*x, sigma_y) 
    ## store data
    tb <- as_data_frame(cbind (x, group, y)) 
    ## fit model
    fit <- lmer(y ~ x + (1 | group))
    se_residual <- get_residual_se(fit)
    sigmaG <- as.data.frame(VarCorr(fit))$sdcor[1]

    tb <- add_pi(tb, fit, includeRanef = includeRanef, alpha = alpha, names = c("lwr_para", "upr_para")) %>%
        add_pi(fit, type = "sim", includeRanef = includeRanef, alpha = alpha, nSims = 1000, names = c("lwr_pi", "upr_pi")) %>%
        add_pi(fit, type = "sim_lme4", includeRanef = includeRanef, alpha = alpha, nSims = 1000, names = c("lwr_lme4", "upr_lme4"))


    y_new_store <- matrix(0, nrow = 10000, ncol = NROW(tb)) 

    for (i in 1:NROW(tb)){

        if (includeRanef)
            y_new_store[,i] <- predict(fit, tb) + rnorm(10000, 0, sd = se_residual)
        else
            y_new_store[,i] <- predict(fit, tb, re.form = NA) + rnorm(10000, 0, sd = se_residual) + rnorm(10000, 0, sd = sigmaG)

        store_para <- store_para + ifelse((tb[i,]$lwr_para < y_new_store) && (y_new_store < tb[i,]$upr_para), 1, 0)
        store_pi <- store_pi + ifelse((tb[i,]$lwr_pi < y_new_store) && (y_new_store < tb[i,]$upr_pi), 1, 0)
        store_para <- store_sim + ifelse((tb[i,]$lwr_sim < y_new_store) && (y_new_store < tb[i,]$upr_sim), 1, 0)
    }
    

    return(list(para_coverage_prob = store_para / sims,
                pi_coverage_prob = store_para / sims,
                boot_coverage_prob = store_para / sims))
}




system.time(results <- sim_study_pi(n = 250, fixed_ef = c(1,1), sigma_y = 1,
                         sigma_g = 1, group_no = 5,
                         sims = 100, bootSims = 200,
                         includeRanef = TRUE, alpha = 0.5)); results
