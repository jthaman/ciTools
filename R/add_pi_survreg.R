# Copyright (C) 2017 Institute for Defense Analyses
#
# This file is part of ciTools.
#
# ciTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ciTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ciTools. If not, see <http://www.gnu.org/licenses/>.

#' Prediction Intervals for Survival Models
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{survreg}.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{survreg}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param method A string. Determines the method used to calculate
#'     prediction intevals. Must be one of either \code{"naive"} or
#'     \code{"boot"}.
#' @param nSims A positive integer. Determines the number of bootstrap
#'     replicates.
#' @param ... Additional arguments.
#'
#' @return A tibble, \code{tb}, with predicted medians, upper and lower
#'     prediction bounds attached.
#'
#' @seealso \code{\link{add_ci.survreg}} for confidence intervals for
#'     \code{survreg} objects, \code{\link{add_probs.survreg}} for
#'     conditional probabilities of \code{survreg} objects, and
#'     \code{\link{add_quantile.survreg}} for response quantiles of
#'     \code{survreg} objects.
#'
#' @examples
#'
#' @export

add_pi.survreg <- function(tb, fit, alpha = 0.05,
                           names = NULL, yhatName = "median_pred",
                           nSims = 10000,
                           method = "naive",
                           ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb)))
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    if (method == "calibrated")
        stop("not yet implemented")
    ## sim_pi_survreg_calibrated(tb, fit, alpha, names, yhatName, nSims)
    else if (method == "naive")
        pi_survreg_naive(tb, fit, alpha, names, yhatName)
    else if (method == "boot")
        sim_pi_survreg_boot(tb, fit, alpha, names, yhatName, nSims)
    else
        stop("unrecognized method.")
}

pi_survreg_naive <- function(tb, fit, alpha, names, yhatName){
    med <- predict(fit, tb, p = 0.5, type = "quantile")
    lwr <- predict(fit, tb, p = alpha / 2, type = "quantile")
    upr <- predict(fit, tb, p = 1 - alpha / 2, type = "quantile")

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- med

    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)
}

## TODO: add argument to return the sim_response matrix for further study?

qsev <- function(p) {
    p=ifelse(p>=.99999999999999,.99999999999999,p)
    p=ifelse(p<=1-.99999999999999,1-.99999999999999,p)
    log(-log(1-p))
}

rsev <- function(n) {
    qsev(runif(n))
}


psev <- function(z) {
    1-exp(-exp(z))
}

dsev <- function(z) {
    exp(z-exp(z))
}

sim_surv_coefs <- function(tb, fit, nSims){
    vcov.hat <- vcov(fit)

    if (fit$dist == "exponential")
        cov_mat <- cbind(rbind(cov_mat, 0), 0)

    beta.logsigma.hat <- c(coef(fit), fit$scale)
    params <- matrix(NA,
                     nrow = nSims,
                     ncol = length(beta.logsigma.hat))

    params <- MASS::mvrnorm(nSims, beta.logsigma.hat, vcov.hat)
    params
}

## only works with log-linear models
get_sim_response_surv_boot <- function(tb, fit, params){
    nSims <- dim(params)[1]
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit, data = tb)
    distr <- fit$dist
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nSims){
        linear_pred <- modmat %*% params[i,1:dim(params)[2] - 1]
        scale <- fit$scale

        if (distr == "lognormal"){
            sim_response[,i] <- exp(linear_pred + scale * rnorm(n = nPreds))
        }
        if (distr == "weibull"){
            sim_response[,i] <- exp(linear_pred + scale * rsev(n = nPreds))
        }
        if (distr == "exponential"){
            sim_response[,i] <- exp(linear_pred + rsev(n = nPreds))
        }
        if (distr == "loglogistic"){
            sim_response[,i] <- exp(linear_pred + scale * rlogis(n = nPreds))
        }
    }
    sim_response
}

sim_pi_survreg_boot <- function(tb, fit, alpha, names, yhatName, nSims){

    out <- predict(fit, newdata = tb, type = "response")

    params <- sim_surv_coefs(tb = tb,
                             fit = fit,
                             nSims = nSims)

    sim_response <- get_sim_response_surv_boot(tb, fit, params)

    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha/2, type = 1)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out

    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)
}

## ### calibrated method ###

## will only work with log-linear models
## currently only works with log normal model
## need to implement a censoring mechanism
## get_sim_response_surv <- function(tb, fit, nSims){
##     nPred <- NROW(tb)
##     beta.hat <- coef(fit)
##     sigma <- fit$scale
##     modmat <- model.matrix(fit, data = tb)
##     dist <- fit$dist
##     sim_response <- matrix(0, ncol = nPred, nrow = nSims)

##     for (i in 1:nSims){
##         linear_pred <- modmat %*% beta.hat
##         if(dist == "lognormal"){
##             sim_response[i,] <- exp(linear_pred + sigma * rnorm(n = nPred))
##         }
##     }
##     sim_response
## }

## get_mle_surv <- function(sim_response, tb, fit){
##     nSims <- NROW(sim_response)
##     dist <- fit$dist
##     for (i in 1:nSims){
##         ## how to get this into data arg? Use a regex...
##         T <- sim_response[i,]
##         form <- formula(fit)
##         temp_model <- survreg(form, data = tb, dist = fit$dist)
##     }
## }

## si_pi_survreg_calibrated <- function(tb, fit, alpha, names, yhatName, nSims){
##     alpha_0 <- alpha
##     sim_dat <- get_sim_response_surv(tb = tb, fit = fit, nSims = nSims)


## }

### boot method ###
