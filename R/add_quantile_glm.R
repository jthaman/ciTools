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
#' compute the quantiles of the response of a Poisson regression with
#' the \eqn{\log}-link function.
#'
#' Quantiles of generalized linear models are determined by
#' \code{add_quantile} through a simulation using \code{arm::sim}.
#'
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
#' @seealso \code{\link{add_ci.lm}} for confidence intervals for
#'     \code{lm} objects, \code{\link{add_pi.lm}} for prediction
#'     intervals of \code{lm} objects, and \code{\link{add_probs.lm}}
#'     for response probabilities of \code{lm} objects.
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
#' @export

add_quantile.glm <- function(tb, fit, p, name = NULL, yhatName = "pred",
                             nSims = 200, ...){
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name))
        name <- paste("quantile", p, sep="")
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }
    if (fit$family$family == "binomial"){
       stop ("Quantiles for Logistic Regression don't make sense") 
    }
    if (fit$family$family == "poisson"){
        warning ("The response is not continuous, so estimated quantiles are only approximate")
        sim_quantile_pois(tb, fit, p, name, yhatName, nSims)
    }
}

sim_quantile_pois <- function(tb, fit, p, name, yhatName, nSims){
    nPreds <- NROW(tb)
    modmat <- model.matrix(fit)
    response_distr <- fit$family$family
    inverselink <- fit$family$linkinv
    out <- predict(fit, newdata = tb, type = "response")
    sims <- arm::sim(fit, n.sims = nSims)
    sim_response <- matrix(0, ncol = nSims, nrow = nPreds)

    for (i in 1:nPreds){
        if(response_distr == "poisson"){
            sim_response[i,] <- rpois(n = nSims, lambda = inverselink(rnorm(nPreds,sims@coef[i,] %*% modmat[i,], sd = sims@sigma[i])))
            }
    }

    quants <- apply(sim_response, 1, FUN = quantile, probs = p, type = 1)

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quants
    tibble::as_data_frame(tb)


}
