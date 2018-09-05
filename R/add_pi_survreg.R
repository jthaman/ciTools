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
#' @param nSims A positive integer. Determines the number of
#' @param bootData A logical. If \code{TRUE} a matrix is bootstrap
#'     replicates is returned.  simulations to run.
#' @param ... Additional arguments.
#'
#' @return A tibble, \code{tb}, with predicted values, upper and lower
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
                           names = NULL, yhatName = "pred",
                           nSims = 2000, bootData = FALSE,
                           ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb)))
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    sim_pi_survreg(tb, fit, alpha, names, yhatName, nSims, bootData)
}

sim_surv_coefs <- function(tb, fit, nSims){
    V.hat <- vcov(fit)
    beta.sigma.hat <- c(coef(fit), fit$scale)
    params <- matrix(NA,
                     nrow = nSims,
                     ncol = length(beta.sigma.hat))

    params <- MASS::mvrnorm(nSims, beta.sigma.hat, V.hat)
    params
}

## TODO
## only work with log-linear models
get_sim_response_surv <- function(tb, fit, params){
    nSims <- dim(params)[1]
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit, data = tb)
    response_distr <- fit$dist
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nSims){
        pred <- modmat %*% params[i,1:dim(params)[2] - 1]
        if(response_distr == "lognormal"){
            sim_response[,i] <- exp(pred  +
                                    params[i, dim(params)[2]] * rnorm(n = nPreds))
        }
    }
    sim_response
}

sim_pi_survreg <- function(tb, fit, alpha, names, yhatName, nSims, bootData){

    out <- predict(fit, newdata = tb, type = "response")

    params <- sim_surv_coefs(tb = tb,
                             fit = fit,
                             nSims = nSims)

    sim_response <- get_sim_response_surv(tb, fit, params)

    if(bootData)
        sim_response
    else {
        lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha/2, type = 1)
        upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)

        if(is.null(tb[[yhatName]]))
            tb[[yhatName]] <- out
        tb[[names[1]]] <- lwr
        tb[[names[2]]] <- upr
        tibble::as_data_frame(tb)
    }
}
