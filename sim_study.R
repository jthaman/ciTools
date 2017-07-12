library(lme4)
library(addr)
library(tidyverse)

## simulation study comparison of LMM CI and PI techniques in addr.

## balanced design data generation

## n = 100;
## fixed_ef = c(1,1);
## sigma_y = 1;
## sigma_g = 1;
## group_no = 5
## sims = 100;
## bootSims = 200;
## includeRanef = TRUE;
## alpha = 0.5

## 100 simulations
system.time(sim_study_ci(n = 100, fixed_ef = c(1,1), sigma_y = 1,
                         sigma_g = 1, group_no = 5,
                         sims = 10, bootSims = 200,
                         includeRanef = TRUE, alpha = 0.5)
            )

sim_study_ci <- function(n, fixed_ef, sigma_y, sigma_g, group_no, sims,
                         bootSims, includeRanef, alpha){

    success_para <- success_pi <- success_boot <- rep(0, sims)

    for (i in 1:sims){
        group <- rep (1:group_no, each = n/group_no) # group labels
        a <- rnorm (group_no, 0, sigma_g) 
        x <- seq(1,2, length.out = n)
        y <- rnorm (n, fixed_ef[1] + a[group] + fixed_ef[2]*x, sigma_y) 
        tb <- as_data_frame(cbind (x, group, y))
        b0 <- fixed_ef[1]
        b1 <- fixed_ef[2]
        EYx <- b0 + b1*x + a[group]
        EYxg <- b0 + b1*x

        if (includeRanef)
            truth <- EYx
        else
            truth <- EYxg

        fit <- lmer(y ~ x + (1 | group))

        tb <- add_ci(tb, fit, includeRanef = includeRanef, alpha = alpha, names = c("lwr_para", "upr_para"))
        tb <- add_ci(tb, fit, type = "sim", includeRanef = TRUE, nSims = 1000, names = c("lwr_pi", "upr_pi"))
        tb <- add_ci(tb,fit, type = "boot", includeRanef = TRUE, nSims = bootSims, names = c("lwr_boot", "upr_boot"))

        if (all( (tb$lwr_para < truth) && (tb$upr_para > truth) ))
            success_para[i] <- 1 
        if (all( (tb$lwr_pi < truth) && (tb$upr_pi > truth) ))
            success_pi[i] <- 1 
        if (all( (tb$lwr_boot < truth) && (tb$upr_boot > truth) ))
            success_boot[i] <- 1 

    }
    return(list(para_coverage_prob = mean(success_para),
                pi_coverage_prob = mean(success_pi),
                boot_coverage_prob = mean(success_boot)))
}















