#' ---
#' title: GLM CI candidates for add_ci
#' author: John Haman
#' date: June 8, 2017
#' output:
#'    html_document:
#'      toc: true
#' ---

#' ## Methods for adding CI to tibble + GLM.
#' choose either parametric (wald) intervals or bootstrap, which
#' should probably not be implemented in the final version. This
#' is basically the same as the method you already made for lmer
#' objects. The default type of interval is "response" meaning the
#' intervals are generated on the response scale. I'm wondering what to
#' call the intervals that can be made at the linear
#' level e.g. "Link level?"

rm(list = ls())
library(tidyverse) 

add_ci <- function(tb, fit, alpha = 0.05, ciNames = c("LCB", "UCB"), ...){
  UseMethod("add_ci")
}

add_ci.glm <- function(tb, fit, alpha = 0.05, ciNames = c("LCB", "UCB"),
                       type = "response", ciType = "parametric", bootReplicates = 1000, ...){
  if(ciType == "bootstrap") return(
    bootstrap_ci_glm(tb, fit, alpha, ciNames, type, bootReplicates, ...)
  ) else
    if(ciType == "parametric") return(
      parametric_ci_glm(tb, fit, alpha, ciNames, type, ...)
    ) else
      if(!(ciType %in% c("bootstrap", "parametric"))) stop("Incorrect interval type specified!")
}

#' Define function to handle parametric intervals. 

parametric_ci_glm <- function(tb, fit, alpha, ciNames, type, ...){
    ## get the fitted values, and the SEs of the fitted
    ## values. Predict for GLMs doesn't give interval estimates
    ## directly. These values are at the link level.
    out <- predict(fit, tb, se.fit = TRUE)
    ## *I have not double checked the degrees of freedom here.*
    crit_val <- qt(p = 1 - alpha/2, df = fit$df.residual)
    ## convenience function
    inverselink <- fit$family$linkinv
    ## get the reponse level intervals
    if(type == "response"){
        upr <- inverselink(out$fit + crit_val * out$se.fit)
        lwr <- inverselink(out$fit - crit_val * out$se.fit)
        pred <- inverselink(out$fit)
    ## get the link level intervals
    }else{
        upr <- out$fit + crit_val * out$se.fit
        lwr <- out$fit - crit_val * out$se.fit
        pred <- out$fit
    }
    ## bind to the tibble and return
    tb[[ciNames[1]]] <- lwr
    if(is.null(tb[["pred"]])) tb[["pred"]] <- pred
    tb[[ciNames[2]]] <- upr
    tb

}

#' Some testing to make sure the generic works okay

parametric_ci_glm(tb, lmod, 0.05, c("LCB", "UCB"), type = "response")
add_ci.glm(tb, lmod, 0.05, c("LCB", "UCB"), type = "response")
add_ci(tb, lmod, 0.05, c("LCB", "UCB"), type = "response")

#' Bootstrap confidence intervals from GLMs
#' need to require(boot) 

## require(boot)

## #' This function is the statistic that the boot() function will apply to resampled data.

## boot_prediction_glm <- function(tb, indices, fit){
##     boot_fit <- glm(formula = fit$formula, data = tb[indices, ], family = fit$family$family)
##     predict(boot_fit, tb)
## }

## #' Or should it be predict(boot_fit, tb[indices,])?
## #' <br>
## #' function to calculate bias corrected confidence intervals, which
## #' may be a better idea than just using the usual non parametric
## #' quantile intervals.  

## get_bca_column <- function(j){
##     boot.ci(boot_store, type = 'bca', index = j)$bca[4:5]
## }

## #' Function to get the boot strap CIs from a GLM 

## bootstrap_ci_glm <- function(tb, fit, alpha, ciNames, type, bootReplicates , ...){
##     ## Perform the bootstrap here
##     boot_store <- boot(data = tb, statistic = boot_prediction_glm, sim = "ordinary", R = bootReplicates, fit = fit)
##     ## This is commented out because BCa intervals are not implemented right now
##     ## ints <- sapply(1:NCOL(boot_store$t), get_bca_column)
##     ## ints <- t(ints)
##     ## upr <- ints[,2]
##     ## lwr <- ints[,1]
##     upr <- apply(boot_store$t, MARGIN = 2, FUN = quantile, probs = 1 - alpha / 2)
##     lwr <- apply(boot_store$t, MARGIN = 2, FUN = quantile, probs = alpha / 2)
##     fitted <- apply(boot_store$t, MARGIN = 2, FUN = quantile, probs = 0.5)
##     inverselink <- fit$family$linkinv
##     if(type == "response"){
##         upr <- inverselink(upr)
##         lwr <- inverselink(lwr)
##         pred <- inverselink(fitted)
##     }else{
##         upr <- upr
##         lwr <- lwr
##         pred <- fitted
##     }
##     tb[[ciNames[1]]] <- lwr
##     if(is.null(tb[["pred"]])) tb[["pred"]] <- pred
##     tb[[ciNames[2]]] <- upr
##     tb
## }

#' ## Testing 

#' Just use the cars data set
tb <- as_data_frame(cars)
#' Poisson Model
pmod <- glm(dist ~ speed, family = poisson, data = cars) 
#' Define some binomial variable from logistic regression
win <- ifelse(cars$dist > 50, 1, 0)
#' Logistic Model
lmod <- glm(win ~ speed, family = binomial, data = tb)

#' Parametric and Bootstrap intervals for Poisson Model
parametric_ci_glm(tb, pmod, 0.05, c("LCB", "UCB"), type = "response")
bootstrap_ci_glm(tb, pmod, bootReplicates = 1000, 0.05, c("LCB", "UCB"), type = "response")

#' Notice that the BS intervals are much bigger... Bigger should be
#' expected I think, because we are leaning less on the assumptions.

#' Parametric and Bootstrap intervals for Logistic Model
parametric_ci_glm(tb, lmod, 0.05, c("LCB", "UCB"), type = "response")
add_ci(tb =tb, fit =lmod, alpha =0.05, ciNames = c("LCB", "UCB"), type = "response")
add_ci.glm(tb, lmod, 0.05, c("LCB", "UCB"), type = "response")

#' Looks reasonable
bootstrap_ci_glm(tb, lmod, bootReplicates = 1000, 0.05, c("LCB", "UCB"), type = "response")

#' This is basically all wrong... but it's not clear to me what the failure is.
#' <br>
#' I'm to open any comments for cleaning up this code.
# rmarkdown::render("c:/Users/jhaman/Documents/R code/add_ci_glm.R")

