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

#' Prediction Intervals for Generalized Linear Models
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{glm}. Prediction Intervals are generated through
#' simulation. At the moment, only prediction intervals for Poisson
#' GLMs with the log link function are supported. Do not try to use
#' this function with a logistic regression, as those prediction
#' intervals are nonsensical.
#'
#' @details Note that if the response is count data, e.g. the response
#'     is modeled with a Poisson regression, prediction intervals are
#'     only approximate.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, prediction bounds will automatically be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param type A string. Currently must be \code{"sim"}.
#' @param nSims A positive integer. Determines the number of
#'     simulations to run.
#' @param ... Additional arguments.
#' 
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' 
#' @export


add_pi.glm <- function(tb, fit, alpha = 0.05, names = NULL,
                       type= "sim", nSims = 1000, ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    }

    if(fit$family$family == "binomial"){
        stop("Prediction interval for Bernoulli response doesn't make sense")
    }

    if(fit$family$family == "poisson"){
        warning("The response is not continuous, so Prediction Intervals are only approximate")
    }

    if(type == "sim"){
        sim_pi_glm(tb, fit, alpha, names, nSims)
    }
    else if(!(type %in% c("sim")))
        stop("Only Simulated prediction intervals are implemented for glm objects")
}

## TODO : hardcode more response distributions
## TODO : Smooth the prediction intervals

sim_pi_glm <- function(tb, fit, alpha, names, nSims){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- inverselink(predict(fit, tb))
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nPreds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims, lambda = inverselink(rnorm(nPreds,sims@coef[i,] %*% modmat[i,], sd = sims@sigma[i])))
            }
    }

    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha / 2, type = 1)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)

    
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)

}

