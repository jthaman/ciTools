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
#' coefficients. At the moment, only prediction intervals for Poisson,
#' Quasipoisson, Gaussian, and Gamma GLMs are supported. Note that if
#' the response is count data, prediction intervals are only
#' approximate. Simulation from the QuasiPoisson model is performed
#' with the negative binomial distribution, see Gelman and Hill
#' (2007).
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param nSims A positive integer. Determines the number of
#'     simulations to run.
#' @param ... Additional arguments.
#'
#' @return A dataframe, \code{df}, with predicted values, upper and lower
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


add_pi.glm <- function(df, fit, alpha = 0.05, names = NULL, yhatName = "pred",
                       nSims = 2000,  ...){

    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(df)))
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")


    if(fit$family$family == "binomial")
      if(max(fit$prior.weights) == 1)
          stop("Prediction intervals for Bernoulli response variables aren't useful")
      else {
          warning("Treating weights as indicating the number of trials for a binomial regression where the response is the proportion of successes")
          warning("The response variable is not continuous so Prediction Intervals are approximate")
        }

    if(fit$family$family %in% c("poisson", "quasipoisson"))
        warning("The response is not continuous, so Prediction Intervals are approximate")

    if(!(fit$family$family %in% c("poisson", "quasipoisson", "Gamma", "binomial", "gaussian")))
        stop("Unsupported family")

    if(fit$family$family == "gaussian")
        pi_gaussian(df, fit, alpha, names, yhatName)
    else
        sim_pi_glm(df, fit, alpha, names, yhatName, nSims)
}

pi_gaussian <- function(df, fit, alpha, names, yhatName){
    sigma_sq <- summary(fit)$dispersion
    inverselink <- fit$family$linkinv
    out <- predict(fit, newdata = df, se.fit = TRUE)
    se_terms <- out$se.fit
    t_quant <- qt(p = alpha/2, df = fit$df.residual, lower.tail = FALSE)
    se_global <- sqrt(sigma_sq + se_terms^2)
    lwr <- inverselink(out$fit) - t_quant * se_global
    upr <- inverselink(out$fit) + t_quant * se_global

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- inverselink(out$fit)
    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}

get_sim_response <- function(df, fit, nSims){

    nPreds <- NROW(df)
    modmat <- model.matrix(fit, data = df)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    overdisp <- summary(fit)$dispersion
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(NA, ncol = nSims, nrow = nPreds)

    for (i in 1:nSims){
        yhat <- inverselink(modmat %*% sims@coef[i,])
        if(response_distr == "poisson"){
            sim_response[,i] <- rpois(n = nPreds,
                                      lambda = yhat)
        }
        if(response_distr == "quasipoisson"){
            a <- yhat / (overdisp - 1)
            sim_response[,i] <- rnegbin(n = nPreds,
                                        mu = yhat,
                                        theta = a)
        }
        if(response_distr == "Gamma"){
            sim_response[,i] <- rgamma(n = nPreds,
                                       shape = 1/overdisp,
                                       rate = 1/(yhat *overdisp))
        }
        if(response_distr == "binomial"){
            yhat <- inverselink(modmat %*% sims@coef[i,]) * fit$prior.weights
            sim_response[,i] <- rbinom(n = nPreds,
                                       size = fit$prior.weights,
                                       prob = yhat / fit$prior.weights)
        }
        if(response_distr == "gaussian"){
            yhat <- inverselink(modmat %*% sims@coef[i,])
            sim_response[,i] <- rnorm(n = nPreds,
                                      mean = yhat,
                                      sd = sqrt(overdisp))
        }

    }
    sim_response
}

sim_pi_glm <- function(df, fit, alpha, names, yhatName, nSims){
    out <- predict(fit, newdata = df, type = "response")
    sim_response <- get_sim_response(df = df, fit = fit, nSims = nSims)
    lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha/2, type = 1)
    upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1)

    if(fit$family$family == "binomial"){
      out <- out * fit$prior.weights
      warning("For binomial models, add_pi's column of fitted values reflect E(Y|X) rather than typical default for logistic regression, pHat")
    }

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- out
    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}
