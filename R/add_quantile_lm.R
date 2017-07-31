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
#' This function is one of the methods of \code{add_quantile}. It is
#' automatically called when \code{add_quantile} is called on objects
#' of class \code{lm}.
#'
#' Quantiles for linear models are determined parametrically, by
#' applying a pivotal quantity to the distribution of \eqn{Y|x}.
#' 
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lm}. Predictions are made with
#'     this object.
#' @param p A real number between 0 and 1. Sets the level of the
#'     quantiles.
#' @param name \code{NULL} or a string. If \code{NULL},
#'     quantiles will automatically be named by \code{add_quantile},
#'     otherwise, they will be named \code{name}
#' @param yhatName A string. Name of the vector of predictions.
#' @param log_response logical. If TRUE, quantiles will be generated
#'     for the prediction made with a log-linear model: \eqn{\log(Y) =
#'     X\beta + \epsilon}. These quantiles will be on the scale of the
#'     original response, \eqn{Y}.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values and level -
#'     \emph{p} quantiles attached.
#'
#' @seealso \code{\link{add_ci.lm}} for confidence intervals for
#'     \code{lm} objects. \code{\link{add_pi.lm}} for prediction
#'     intervals of \code{lm} objects, and \code{\link{add_probs.lm}}
#'     for response probabilities of \code{lm} objects.
#'
#' @examples
#'
#' # Fit a linear Model
#' fit <- lm(dist ~ speed, data = cars)
#'
#' # Find the 0.7-quantile (70th percentile) of new dists, given the linear model fit.
#' add_quantile(cars, fit, p = 0.7)
#'
#' # As above, but with a custom name for the vector of quantiles
#' add_quantile(cars, fit, p = 0.7, name = "my_quantile")
#' 
#' @export

add_quantile.lm <- function(tb, fit, p, name = NULL, yhatName = "pred",
                            log_response = FALSE, ...){
    if (log_response)
        add_quantile_lm_log(tb, fit, p, name, yhatName)
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
        if (is.null(tb[[yhatName]]))
            tb[[yhatName]] <- fitted
        tb[[name]] <- out_quantiles
        tibble::as_data_frame(tb)
    }
}

