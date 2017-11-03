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

#' Quantiles for the Response of a Generalized Linear Model
#'
#' This function is one of the methods of
#' \code{add_quantile}. Currently, you can only use this function to
#' compute the quantiles of the response of Poisson, Quasipoisson, or
#' Gamma regression models.  Quantile estimates for Bernoulli response
#' variables (i.e., logistic regression) are not supported.
#' 
#' Quantiles of generalized linear models are determined by
#' \code{add_quantile} through a simulation using \code{arm::sim}. If
#' a Quasipoisson regression model is fit, simulation using the
#' Negative Binomial distribution is performed, see Gelman and Hill
#' (2007).
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{glm}. Predictions are made with
#'     this object.
#' @param p A real number between 0 and 1. Sets the probability level
#'     of the quantiles.
#' @param name \code{NULL} or a string. If \code{NULL},
#'     quantiles automatically will be named by \code{add_quantile},
#'     otherwise, they will be named \code{name}.
#' @param yhatName A string. Name of the vector of predictions.
#' @param nSims A positive integer. Set the number of simulated draws
#'     to use.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values and level
#'     \emph{p} quantiles attached.
#'
#' @seealso \code{\link{add_ci.glm}} for confidence intervals for
#'     \code{glm} objects, \code{\link{add_pi.glm}} for prediction
#'     intervals of \code{glm} objects, and \code{\link{add_probs.glm}}
#'     for response probabilities of \code{glm} objects.
#'
#' @examples
#'
#' # Fit a Poisson GLM
#' fit <- glm(dist ~ speed, data = cars, family = "poisson")
#'
#' # What is the 0.3-quantile (or 30th percentile) of new distances,
#' # given the Poisson model?
#' add_quantile(cars, fit, p = 0.3)
#'
#' # As above, but now find the 0.5-quantile (50th percentile), change
#' # the number of simulations to run, and give the vector of
#' # quantiles a custom name.
#' add_quantile(cars, fit, p = 0.5, name = "my_quantile", nSims = 300)
#'
#' 
#' @export

add_quantile.glm <- function(tb, fit, p, name = NULL, yhatName = "pred",
                             nSims = 2000, ...){
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name))
        name <- paste("quantile", p, sep="")
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }

    if (fit$family$family == "binomial"){
        if(max(fit$prior.weights) == 1)
        stop("Prediction intervals for Bernoulli response variables aren't useful") else {
          warning("Treating weights as indicating the number of trials for a binomial regression where the response is the proportion of successes")
          warning("The response variable is not continuous so Prediction Intervals are approximate")
        }
    }

    if (fit$family$family %in% c("poisson", "qausipoisson"))
        warning("The response is not continuous, so estimated quantiles are only approximate")

    if(!(fit$family$family %in% c("poisson", "quasipoisson", "Gamma", "binomial")))
        stop("Unsupported family")

    sim_quantile_other(tb, fit, p, name, yhatName, nSims)
}

sim_quantile_other <- function(tb, fit, p, name, yhatName, nSims){

    out <- predict(fit, newdata = tb, type = "response")
    sim_response <- get_sim_response(tb, fit, nSims)
    quants <- apply(sim_response, 1, FUN = quantile, probs = p, type = 1)

    if(fit$family$family == "binomial"){
      out <- out * fit$prior.weights
      warning("For binomial models, add_quantile's column of fitted values reflect E(Y|X) rather than typical default for logistic regression, pHat")
    }
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quants
    tibble::as_data_frame(tb)
}
