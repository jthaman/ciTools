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

#' Response Probabilities for Negative Binomial Models
#'
#' This is the method \code{add_probs} uses if the model fit is an
#' object of class \code{negbin}. Probabilities are determined through
#' simulation, using the same method as \code{add_pi.negbin}.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{negbin}. Predictions are made with
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
#' @seealso \code{\link{add_ci.negbin}} for confidence intervals for
#'     \code{negbin} objects, \code{\link{add_pi.negbin}} for prediction
#'     intervals of \code{negbin} objects, and
#'     \code{\link{add_quantile.negbin}} for response quantiles of
#'     \code{negbin} objects.
#'
#' @examples
#' 
#' @export
add_probs.negbin <- function(tb, fit, q, name = NULL, yhatName = "pred",
                             comparison = "<", nSims = 2000, ...){

    if (is.null(name) & (comparison == "<"))
        name <- paste("prob_less_than", q, sep="")
    if (is.null(name) & (comparison == ">"))
        name <- paste("prob_greater_than", q, sep="")
    if (is.null(name) & (comparison == "<="))
        name <- paste("prob_less_than_or_equal_to", q, sep="")
    if (is.null(name) & (comparison == ">="))
        name <- paste("prob_greater_than_or_equal_to", q, sep="")
    if (is.null(name) & (comparison == "="))
        name <- paste("prob_equal_to", q, sep="")

    if ((name %in% colnames(tb)))
        warning ("These probabilities may have already been appended to your dataframe. Overwriting.")

    sim_probs_negbin(tb, fit, q, name, yhatName, nSims, comparison)
}


sim_probs_negbin <- function(tb, fit, q, name, yhatName, nSims, comparison){

    out <- predict(fit, newdata = tb, type = "response")
    sim_response <- get_sim_response_nb(tb, fit, nSims)

    probs <- apply(sim_response, 1, FUN = calc_prob, quant = q, comparison = comparison)
    
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- probs
    tibble::as_data_frame(tb)
}

