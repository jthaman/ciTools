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
#' class \code{glm}.
#'
#' Prediction intervals are generated through simulation with the aid
#' \code{arm::sim}, which simulates the uncertainty in the regression
#' coefficients. At the moment, only prediction intervals for Poisson
#' and Quasipoisson GLMs is supported. Note that if the response is
#' count data, prediction intervals are only approximate.
#' 
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param type A string. Currently \code{type = "sim"} is the only
#'     valid string.
#' @param nSims A positive integer. Determines the number of
#'     simulations to run.
#' @param ... Additional arguments.
#' 
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @seealso \code{\link{add_ci.glm}} for confidence intervals for
#'     \code{glm} objects, \code{\link{add_probs.glm}} for conditional
#'     probabilities of \code{glm} objects, and
#'     \code{\link{add_quantile.glm}} for response quantiles of
#'     \code{glm} objects.
#'
#' @examples
#' # Fit a Poisson model
#' fit <- glm(dist ~ speed, data = cars, family = "poisson")
#' # Add prediction intervals and fitted values to the original data frame
#' add_pi(cars, fit)
#' # Try a different confidence level
#' add_pi(cars, fit, alpha = 0.5)
#' # Try custom names for the prediction bounds (may be useful for plotting)
#' add_pi(cars, fit, alpha = 0.5, names = c("lwr", "upr"))
#' 
#' @export


add_pi.glm <- function(tb, fit, alpha = 0.05, names = NULL, yhatName = "pred", 
                       nSims = 2000, type = "sim", ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) 
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    
    if(fit$family$family == "binomial")
        stop("Prediction interval for Bernoulli response doesn't make sense")
    
    if(fit$family$family %in% c("poisson", "quasipoisson"))
        warning("The response is not continuous, so Prediction Intervals are only approximate")

    if(type == "sim")
        sim_pi_glm(tb, fit, alpha, names, yhatName, nSims)
    
    else if(!(type %in% c("sim")))
        stop("Only Simulated prediction intervals are implemented for glm objects")
}

## We might be able to get a more general solution if we use a
## bootstrap

sim_pi_glm <- function(tb, fit, alpha, names, yhatName, nSims){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- inverselink(predict(fit, tb))
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)
    overdisp <- summary(fit)$dispersion

    for (i in 1:nPreds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims,
                                      lambda = inverselink(sims@coef[i,] %*% modmat[i,]))
        }
        if(response_distr == "quasipoisson"){
            a <- inverselink (modmat[i,] %*% sims@coef[i,]) / (overdisp - 1)
            sim_response[i,] <- rnegbin(n = nSims,
                                        mu = inverselink(sims@coef[i,] %*% modmat[i,]),
                                        theta = a)
        }
        if(response_distr == "Gamma"){
            sim_response[i,] <- rgamma(n = nSims,
                                       shape = 1/overdisp,
                                       rate = 1/inverselink(sims@coef[i,] %*% modmat[i,]) * 1/overdisp)
        }
    }

    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha / 2, type = 1)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)

}

