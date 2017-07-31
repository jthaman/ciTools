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

#' Confidence Intervals for Generalized Linear Model Predictions
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{glm}. Confidence Intervals are determined by making an
#' interval on the scale of the linear predictor, then applying the
#' inverse link function from the model fit to transform the linear
#' level confidence intervals to the response level. 
#'
#' @param tb A tibble or data frame of new data
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the vector of predictions
#'     made for each observation in tb
#' @param type A string, currently \code{type = "parametric"} is the
#'     only supported method.
#' @param response A logical. The default is \code{TRUE}. If
#'     \code{TRUE}, the confidence intervals will be determined for
#'     the expected response, if \code{FALSE}, confidence intervals
#'     will be made on the scale of the linear predictor.
#' @param ... Additional arguments.
#' 
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @seealso \code{\link{add_pi.glm}} for prediction intervals for
#'     \code{glm} objects. \code{\link{add_probs.glm}} for conditional
#'     probabilities of \code{glm} objects, and
#'     \code{\link{add_quantile.glm}} for response quantiles of
#'     \code{glm} objects.
#'
#' @examples
#' # Poisson Regression
#' fit <- glm(dist ~ speed, data = cars, family = "poisson")
#' add_ci(cars, fit)
#' # Try a different confidence level
#' add_ci(cars, fit, alpha = 0.5)
#' # Add custom names to the confidence bounds (may be useful for plotting)
#' add_ci(cars, fit, alpha = 0.5, names = c("lwr", "upr"))
#' 
#' # Logistic Regression
#' fit2 <- glm(I(dist > 30) ~ speed, data = cars, family = "binomial")
#' dat <- cbind(cars, I(cars$dist > 30))
#' add_ci(dat, fit)
#' add_ci(dat, fit, alpha = 0.5)
#' # Make confidence intervals on the scale of the linear predictor
#' add_ci(dat, fit, alpha = 0.5, response = FALSE)
#' # Add custom names to the confidence bounds
#' add_ci(dat, fit, alpha = 0.5, names = c("lwr", "upr"))
#'
#'
#' @export

add_ci.glm <- function(tb, fit, alpha = 0.05, names = NULL, yhatName = "pred",
                      response = TRUE, type = "parametric", ...){

    if (!(fit$converged))
        warning ("coverage probabilities may be inaccurate if the model did not converge")
    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe. Overwriting.")
    }
    if (type == "boot")
        stop ("not yet implemented")
    else if (type == "parametric")
        parametric_ci_glm(tb, fit, alpha, names, yhatName, response)
    else
        if(!(type %in% c("boot", "parametric"))) stop("Incorrect interval type specified!")

}

parametric_ci_glm <- function(tb, fit, alpha, names, yhatName, response){
    out <- predict(fit, tb, se.fit = TRUE)

    crit_val <- qt(p = 1 - alpha/2, df = fit$df.residual)
    inverselink <- fit$family$linkinv
    if (response){
        upr <- inverselink(out$fit + crit_val * out$se.fit)
        lwr <- inverselink(out$fit - crit_val * out$se.fit)
        pred <- inverselink(out$fit)
    }else{
        upr <- out$fit + crit_val * out$se.fit
        lwr <- out$fit - crit_val * out$se.fit
        pred <- out$fit
    }
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- pred
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)

}
