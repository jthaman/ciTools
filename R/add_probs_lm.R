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

#' Response Level Probabilities for Linear Models
#'
#' This is the method \code{add_probs} uses if the model is of class
#' \code{lm}. Probabilities are calculated parametrically,
#' using a pivotal quantity.
#'
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{lm}. Predictions are made with
#'     this object.
#' @param q A real number. A quantile of the response distribution.
#' @param name \code{NULL} or a string. If \code{NULL}, probabilities
#'     automatically will be named by \code{add_probs}, otherwise, the
#'     probabilities will be named \code{name} in the returned data frame.
#' @param yhatName A character vector of length one. Names of the
#' @param comparison \code{"<"}, or \code{">"}. If \code{comparison =
#'     "<"}, then \eqn{Pr(Y|x < q)} is calculated for each observation in
#'     \code{df}. Otherwise, \eqn{Pr(Y|x > q)} is calculated.
#' @param log_response A logical. Default is \code{FALSE}. Set to
#'     \code{TRUE} if the model is log-linear: \eqn{\log(Y) = X \beta
#'     + \epsilon}.
#' @param ... Additional arguments.
#'
#' @return A dataframe, \code{df}, with predicted values and
#'     probabilities attached.
#'
#' @seealso \code{\link{add_ci.lm}} for confidence intervals for
#'     \code{lm} objects, \code{\link{add_pi.lm}} for prediction
#'     intervals of \code{lm} objects, and
#'     \code{\link{add_quantile.lm}} for response quantiles of
#'     \code{lm} objects.
#'
#' @examples
#'
#' # Fit a linear model
#' fit <- lm(dist ~ speed, data = cars)
#'
#' # Calculate the probability that a new dist will be less than 20,
#' # given the model.
#' add_probs(cars, fit, q = 20)
#'
#' # Calculate the probability that a new dist will be greater than
#' # 30, given the model.
#' add_probs(cars, fit, q = 30, comparison = ">")
#'
#' @export

add_probs.lm <- function(df, fit, q, name = NULL, yhatName = "pred",
                         comparison = "<", log_response = FALSE, ...){

    if (is.null(name) && comparison == "<")
        name <- paste("prob_less_than", q, sep="")
    if (is.null(name) && comparison == ">")
        name <- paste("prob_greater_than", q, sep="")
    if ((name %in% colnames(df))) {
        warning ("These probabilities may have already been appended to your dataframe. Overwriting.")
    }

    if (log_response)
        q <- log(q)

    out <- predict(fit, df, interval = "prediction", se.fit = TRUE)
    fitted <- out$fit[,1]
    residual_df <- out$df
    se_fitted <- out$se.fit
    resid_var <- out$residual.scale^2
    se_pred <- sqrt(resid_var + se_fitted^2)
    t_quantile <- (q - fitted) / se_pred
    if (comparison == "<")
        t_prob <- pt(q = t_quantile, df = residual_df)
    if (comparison == ">")
        t_prob <- 1 - pt(q = t_quantile, df = residual_df)
    if (is.null(df[[yhatName]]))
        df[[yhatName]] <- fitted
    df[[name]] <- t_prob
    data.frame(df)
}
