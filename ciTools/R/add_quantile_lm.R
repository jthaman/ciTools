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

#' Quantiles for the Response of a Linear Model
#'
#' This function is one of the methods for \code{add_quantile}. It is
#' automatically called when \code{add_quantile} is applied to objects
#' of class \code{lm}.
#'
#' Quantiles for linear models are determined parametrically, more or
#' less the same way that prediction intervals for this class of
#' models is determined.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param p A real number between 0 and 1. Sets the probability for
#'     which to calculate the quantiles.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     quantiles will automatically be named by \code{add_quantile},
#'     otherwise, they will be named \code{name}
#' @param log_response logical. If TRUE, quantiles will be generated
#'     for the prediction made with a log-linear model: \eqn{\log(Y) =
#'     X\beta + \epsilon}. These quantiles will be on the scale of the
#'     original response, Y.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values and level -
#'     \emph{p} quantile attached.
#'
#' @seealso \code{{\link{add_ci.lm}}} for confidence intervals for
#'     \code{lm} objects. \code{\link{add_pi.lm}} for prediction
#'     intervals of \code{lm} objects, and \code{\link{add_probs.lm}}
#'     for response probabilities of \code{lm} objects.
#'
#' @examples
#' fit <- lm(dist ~ speed, data = cars)
#' add_quantile(cars, fit, p = 0.7)
#' add_quantile(cars, fit, p = 0.7, name = "my_quantile")
#' 
#' @export

add_quantile.lm <- function(tb, fit, p, name = NULL, log_response = FALSE, ...){
    if (log_response)
        add_quantile_lm_log(tb, fit, p, name)
    else {
        if (p == 0.5)
            warning ("The 0.5 quantile is equal to the fitted values")
        if (p <= 0 || p >= 1)
            stop ("p should be in (0,1)")
        if (is.null(name))
            name <- paste("quantile", p, sep="")
        if (name %in% colnames(tb)) {
            warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
        }
        out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
        fitted <- out$fit[,1]
        residual_df <- out$df
        se_fitted <- out$se.fit
        resid_var <- out$residual.scale^2
        se_pred <- sqrt(resid_var + se_fitted^2)
        t_quantile <- qt(p = p, df = residual_df)
        out_quantiles <- fitted + se_pred * t_quantile
        if (is.null(tb[["pred"]]))
            tb[["pred"]] <- fitted
        tb[[name]] <- out_quantiles
        tibble::as_data_frame(tb)
    }
}

