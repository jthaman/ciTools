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

#' Prediction Intervals for Accelerated Failure Time Models
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{survreg}.
#'
#' \code{add_pi.survreg} creates prediction intervals for the survival
#' time $T$ conditioned on the covariates of the \code{survreg}
#' model. In simple terms, this function calculates error bounds
#' within which one can expect to observe a new survival time. Like
#' other parametric survival methods in \code{ciTools}, prediction
#' intervals are limited to unweighted lognormal, exponential,
#' weibull, and loglogistic AFT models.
#'
#' Two methods are available for creating prediction intervals, the
#' "naive" method (Meeker and Escobar, chapter 8) and a simulation
#' method that implements a parametric bootstrap routine. The "naive"
#' method calculates quantiles of the fitted survival time
#' distribution to determine prediction intervals. The parametric
#' bootstrap method simulates new survival times from the conditional
#' survival time distribution, taking into account the uncertainty in
#' the regression coefficients. The bootstrap method is similar to the
#' one implemented in \code{add_pi.glm}.
#'
#' Note: Due to a limitation, the \code{Surv} object must be specified in
#' \code{survreg} function call. See the examples section for one way
#' to do this.
#'
#' Note: \code{add_pi.survreg} cannot inspect the convergence of
#' \code{fit}. Poor maximum likelihood estimates will result in poor
#' prediction intervals. Inspect any warning messages given from
#' \code{survreg}.
#'
#' @param df A data frame of new data.
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
#'     prediction intervals. Must be one of either \code{"naive"} or
#'     \code{"boot"}.
#' @param nSims A positive integer. Determines the number of bootstrap
#'     replicates if \code{method = "boot"}.
#' @param ... Additional arguments.
#'
#' @return A dataframe, \code{df}, with predicted medians, upper and lower
#'     prediction bounds attached.
#'
#' @seealso \code{\link{add_ci.survreg}} for confidence intervals for
#'     \code{survreg} objects, \code{\link{add_probs.survreg}} for
#'     conditional survival probabilities of \code{survreg} objects, and
#'     \code{\link{add_quantile.survreg}} for survival time quantiles
#'     of \code{survreg} objects.
#'
#' @references
#' For a discussion prediction intervals of accelerated failure time models:
#' Meeker, William Q., and Luis A. Escobar. Statistical methods for reliability data. John Wiley & Sons, 2014. (Chapter 8)
#'
#' @examples
#' ## Define a data set.
#' df <- survival::stanford2
#' ## remove a covariate with missing values.
#' df <- df[, 1:4]
#' ## next, create the Surv object inside the survreg call:
#' fit <- survival::survreg(survival::Surv(time, status) ~ age + I(age^2),
#'                          data = df, dist = "lognormal")
#' add_pi(df, fit, alpha = 0.1, names = c("lwr", "upr"))
#'
#' ## Try a different model:
#' fit2 <- survival::survreg(survival::Surv(time, status) ~ age + I(age^2),
#'                           data = df, dist = "weibull")
#' add_pi(df, fit2, alpha = 0.1, names = c("lwr", "upr"))
#'
#' @export

add_pi.survreg <- function(df, fit, alpha = 0.05,
                           names = NULL,
                           yhatName = "median_pred",
                           nSims = 10000,
                           method = "naive",
                           ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }

    if ((names[1] %in% colnames(df)))
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")

    if (!(fit$dist %in%
          c("loglogistic", "lognormal", "loggaussian", "exponential", "weibull")))
        stop("Unsupported distribution")

    if (!is.null(fit$weights))
        if (var(fit$weights) != 0)
            stop("weighted regression is unsupported.")

    if(any(is.na(df)))
        stop("Check df for missingness")

    if (method == "naive")
        pi_survreg_naive(df, fit, alpha, names, yhatName)
    else if (method == "boot")
        sim_pi_survreg_boot(df, fit, alpha, names, yhatName, nSims)
    else
        stop("unrecognized method.")
}

pi_survreg_naive <- function(df, fit, alpha, names, yhatName){
    med <- predict(fit, df, p = 0.5, type = "quantile")
    lwr <- predict(fit, df, p = alpha / 2, type = "quantile")
    upr <- predict(fit, df, p = 1 - alpha / 2, type = "quantile")

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- med

    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}

## Loglogistic distribution functions (taken from SPREDA library)
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

sim_surv_coefs <- function(df, fit, nSims){
    vcov.hat <- vcov(fit)
    beta.hat <- coef(fit)
    params <- matrix(NA,
                     nrow = nSims,
                     ncol = length(beta.hat))

    params <- MASS::mvrnorm(nSims, beta.hat, vcov.hat)
    params
}

get_sim_response_surv_boot <- function(df, fit, params){
    nSims <- dim(params)[1]
    nPreds <- NROW(df)
    modmat <- model.matrix(fit, data = df)
    distr <- fit$dist
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)
    scale <- fit$scale

    for (i in 1:nSims){
        linear_pred <- modmat %*% params[i,1:dim(params)[2]]

        if ((distr == "lognormal") || (distr == "loggaussian")){
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

sim_pi_survreg_boot <- function(df, fit, alpha, names, yhatName, nSims){

    out <- predict(fit, newdata = df, type = "quantile", p = 0.5)

    params <- sim_surv_coefs(df = df,
                             fit = fit,
                             nSims = nSims)

    sim_response <- get_sim_response_surv_boot(df, fit, params)

    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha/2, type = 1)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- out

    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}
