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

#' Response Probabilities for Generalized Linear Models
#'
#' This is the method \code{add_probs} uses if the model fit is an
#' object of class \code{glm}. Probabilities are determined through
#' simulation, using the same method as \code{add_pi.glm}. Currently,
#' only logistic, Poisson, Quasipoisson, and Gamma models are
#' supported.
#'
#' Any of the five comparisons may be made for a Poisson,
#' quasipoisson, or Gamma model: \code{comparison = "<"}, \code{">"},
#' \code{"="}, \code{"<="}, or \code{">="}. For logistic regression,
#' the comparison statement must be equivalent to \eqn{Pr(Y|x = 0)} or
#' \eqn{Pr(Y|x = 1)}.
#'
#' If \code{add_probs} is called on a Poisson, quasiPoisson or Gamma
#' model, a simulation is performed using \code{arm::sim}.
#'
#' If \code{add_probs} is called on a logistic model, the fitted
#' probabilities are used directly (no simulation is required).
#' 
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{glm}. Predictions are made with
#'     this object.
#' @param q A double. A quantile of the response distribution.
#' @param name \code{NULL} or a string. If \code{NULL}, probabilities
#'     automatically will be named by \code{add_probs()}, otherwise,
#'     the probabilities will be named \code{name} in the returned
#'     tibble
#' @param yhatName A string. Name of the vector of predictions.
#' @param comparison A character vector of length one. If
#'     \code{comparison = "<"}, then \eqn{Pr(Y|X < q)} is
#'     calculated. Any comparison is allowed in Poisson regression,
#'     but only certain comparisons may be made in Logistic
#'     regression. See the Details section.
#' @param nSims A positive integer. Controls the number of simulated
#'     draws to make if the model is Poisson.
#' @param ... Additional arguments.
#' 
#' @return A tibble, \code{tb}, with predicted values and
#'     probabilities attached.
#' 
#' @seealso \code{\link{add_ci.glm}} for confidence intervals for
#'     \code{glm} objects, \code{\link{add_pi.glm}} for prediction
#'     intervals of \code{glm} objects, and
#'     \code{\link{add_quantile.glm}} for response quantiles of
#'     \code{glm} objects.
#'
#' @examples
#' # Fit a Poisson model
#' fit <- glm(dist ~ speed, data = cars, family = "poisson")
#'
#' # Determine the probability that a new dist is less than 20, given
#' # the Poisson model.
#' add_probs(cars, fit, q = 20)
#'
#' # Determine the probability that a new dist is greater than 20,
#' # given the Poisson model.
#' add_probs(cars, fit, q = 30, comparison = ">")
#'
#' # Determine the probability that a new dist is greater than or
#' # equal to 20, given the Poisson model.
#' add_probs(cars, fit, q = 30, comparison = ">=")
#'
#' # Fit a logistic model
#' fit2 <- glm(I(dist > 30) ~ speed, data = cars, family = "binomial")
#' add_probs(cars, fit2, q = 0, comparison = "=")
#' add_probs(cars, fit2, q = 1, comparison = "=")
#' 
#' @export
add_probs.glm <- function(tb, fit, q, name = NULL, yhatName = "pred",
                          comparison = "<", nSims = 2000, ...){
    if (is.null(name) && comparison == "<")
        name <- paste("prob_less_than", q, sep="")
    else if (is.null(name) && comparison == ">")
        name <- paste("prob_greater_than", q, sep="")
    else if (is.null(name) && comparison == "<=")
        name <- paste("prob_less_than_or_equal", q, sep="")
    else if (is.null(name) && comparison == ">=")
        name <- paste("prob_greater_than_or_equal", q, sep="")
    else if (is.null(name) && comparison == "=")
        name <- paste("prob_equal_to", q, sep="")

    if ((name %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe. Overwriting.")
    }

    if (fit$family$family == "binomial"){
      if(max(fit$prior.weights) == 1){
        warning("Equivalent to Pr(Y = 0) (or Pr(Y = 1) if comparison = '>' is specified)")
        probs_logistic(tb, fit, q, name, yhatName, comparison)
      } 
      else {
          warning("Treating weights as indicating the number of trials for a binomial regression where the response is the proportion of successes")
          sim_probs_other(tb, fit, q, name, yhatName, nSims, comparison)
      }
    }
    if (fit$family$family %in% c("poisson", "qausipoisson"))
        warning("The response is not continuous, so estimated probabilities are only approximate")

    if (fit$family$family %in% c("poisson", "qausipoisson", "Gamma"))
        sim_probs_other(tb, fit, q, name, yhatName, nSims, comparison)
}

probs_logistic <- function(tb, fit, q, name, yhatName, comparison){
    inverselink <- fit$family$linkinv
    out <- predict(fit, tb, se.fit = TRUE)
    out <- inverselink(out$fit)
    if (((comparison == "=") && (q == 0)) || ((comparison == "<") && (q < 1) && (q > 0)))
        probs <- 1 - out
    else
        probs <- out
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- probs
    tibble::as_data_frame(tb)
}

sim_probs_other <- function(tb, fit, q, name, yhatName, nSims, comparison){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit, data = tb)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)
    overdisp <- summary(fit)$dispersion

    for (i in 1:nSims){

        yhat <- inverselink(modmat %*% sims@coef[i,])
        if(response_distr == "poisson"){
            sim_response[,i] <- rpois(n = nPreds,
                                      lambda = yhat)
        }
        if(response_distr == "quasipoisson"){
            a <- inverselink (modmat %*% sims@coef[i,]) / (overdisp - 1)
            sim_response[,i] <- rnegbin(n = nPreds,
                                        mu = yhat,
                                        theta = a)
        }
        if(response_distr == "Gamma"){
            sim_response[,i] <- rgamma(n = nPreds,
                                       shape = 1/overdisp,
                                       rate = 1/yhat * 1/overdisp)
        }
        if(response_distr == "binomial"){
            yhat <- yhat * fit$prior.weights 
            sim_response[,i] <- rbinom(n = nPreds, 
                                       size = fit$prior.weights,
                                       p = yhat)
        }
    }

    probs <- apply(sim_response, 1, FUN = calc_prob, quant = q, comparison = comparison)
    
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- probs
    tibble::as_data_frame(tb)

}


