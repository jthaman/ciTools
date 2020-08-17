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

#' Quantiles for the Response of a Negative Binomial Regression
#'
#' This function is one of the methods of \code{add_quantile}.
#' 
#' Quantiles of Negative Binomial linear models are determined by
#' \code{add_quantile} through a simulation using \code{arm::sim}.
#' 
#' @param df A data frame of new data.
#' @param fit An object of class \code{negbin}. Predictions are made with
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
#' @return A dataframe, \code{df}, with predicted values and level
#'     \emph{p} quantiles attached.
#'
#' @seealso \code{\link{add_ci.negbin}} for confidence intervals for
#'     \code{negbin} objects, \code{\link{add_pi.negbin}} for prediction
#'     intervals of \code{negbin} objects, and \code{\link{add_probs.negbin}}
#'     for response probabilities of \code{negbin} objects.
#'
#' @examples
#' x1 <- rnorm(100, mean = 1)
#' y <- MASS::rnegbin(n = 100, mu = exp(1 + x1), theta = 5)
#' df <- data.frame(x1 = x1, y = y)
#' fit <- MASS::glm.nb(y ~ x1, data = df)
#' add_quantile(df, fit, p = 0.3)
#' 
#' @export

add_quantile.negbin <- function(df, fit, p, name = NULL, yhatName = "pred",
                                nSims = 2000, ...){
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name))
        name <- paste("quantile", p, sep="")
    if ((name %in% colnames(df))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }

    sim_quantile_nb(df, fit, p, name, yhatName, nSims)
}

sim_quantile_nb <- function(df, fit, p, name, yhatName, nSims){

    out <- predict(fit, newdata = df, type = "response")
    sim_response <- get_sim_response_nb(df, fit, nSims)
    quants <- apply(sim_response, 1, FUN = quantile, probs = p, type = 1)

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- out
    df[[name]] <- quants
    data.frame(df)
}
