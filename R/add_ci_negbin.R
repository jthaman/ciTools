## Copyright (C) 2017 Institute for Defense Analyses
##
## This file is part of ciTools.
##
## ciTools is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## ciTools is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with ciTools. If not, see <http://www.gnu.org/licenses/>.

#' Confidence Intervals for Negative Binomial Linear Model Predictions
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{negbin}.
#'
#' The default link function is log-link. Confidence Intervals are
#' determined by making an interval on the scale of the linear
#' predictor, then applying the inverse link function from the model
#' fit to transform the linear level confidence intervals to the
#' response level. Alternatively, bootstrap confidence intervals may
#' be formed. The bootstrap intervals are formed by first resampling
#' cases from the data frame used to calculate \code{fit}, then bias
#' corrected and accelerated intervals are calculated. See
#' \code{boot::boot.ci} for more details.
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{negbin}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds automatically will be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the vector of predictions made
#'     for each observation in df
#' @param type A string. Must be \code{type = "parametric"} or
#'     \code{type = "boot"}. \code{type} determines the method used to
#'     compute the confidence intervals.
#' @param response A logical. The default is \code{TRUE}. If
#'     \code{TRUE}, the confidence intervals will be determined for
#'     the expected response; if \code{FALSE}, confidence intervals
#'     will be made on the scale of the linear predictor.
#' @param nSims An integer. Number of simulations to perform if the
#'     bootstrap method is used.
#' @param ... Additional arguments.
#'
#' @return A dataframe, \code{df}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @seealso \code{\link{add_pi.negbin}} for prediction intervals for
#'     \code{negbin} objects, \code{\link{add_probs.negbin}} for conditional
#'     probabilities of \code{negbin} objects, and
#'     \code{\link{add_quantile.negbin}} for response quantiles of
#'     \code{negbin} objects.
#'
#' @examples
#' x1 <- rnorm(300, mean = 1)
#' y <- MASS::rnegbin(n = 300, mu = exp(5 + 0.5 * x1), theta = 2)
#' df <- data.frame(x1 = x1, y = y)
#' fit <- MASS::glm.nb(y ~ x1, data = df)
#' df <- df[sample(100),]
#' add_ci(df, fit, names = c("lcb", "ucb"))
#'
#' @export

add_ci.negbin <- function(df, fit, alpha = 0.05, names = NULL, yhatName = "pred",
                          response = TRUE, type = "parametric", nSims = 2000, ...){

    if (!(fit$converged))
        warning ("coverage probabilities may be inaccurate if the model did not converge")
    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(df))) {
        warning ("These CIs may have already been appended to your dataframe. Overwriting.")
    }

    if (type == "boot")
        boot_ci_negbin(df, fit, alpha, names, yhatName, response, nSims)
    else if (type == "parametric")
        parametric_ci_negbin(df, fit, alpha, names, yhatName, response)
    else
        stop("Incorrect interval type specified!")

}

parametric_ci_negbin <- function(df, fit, alpha, names, yhatName, response){
    out <- predict(fit, df, se.fit = TRUE, type = "link")

    crit_val <- qt(p = 1 - alpha/2, df = fit$df.residual)

    inverselink <- fit$family$linkinv

    if (response){
        pred <- inverselink(out$fit)
        upr <- inverselink(out$fit + crit_val * out$se.fit)
        lwr <- inverselink(out$fit - crit_val * out$se.fit)
        if(fit$family$link %in% c("inverse", "1/mu^2")){
            upr1 <- lwr
            lwr <- upr
            upr <- upr1
        }
    }else{
        upr <- out$fit + crit_val * out$se.fit
        lwr <- out$fit - crit_val * out$se.fit
        pred <- out$fit
    }

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- pred
    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}


boot_fit_nb<- function(data, df, fit, lvl, indices){
    data_temp <- data[indices,]
    form <- fit$call[2]
    ##temp_fit <- glm.nb(form, data = data_temp)
    temp_fit <- update(fit, data = data_temp)
    predict(temp_fit, newdata = df, type = lvl)
}

boot_ci_negbin <- function(df, fit, alpha, names, yhatName, response, nSims){

    if (response){
        lvl <- "response"
    }
    else{
        lvl <- "link"
    }

    out <- predict(fit, df, type = lvl)

    boot_obj <- boot(data = fit$model,
                     statistic = boot_fit_nb,
                     R = nSims,
                     fit = fit,
                     df = df,
                     lvl = lvl)


    temp_mat <- matrix(0, ncol = 2, nrow = NROW(df))

    for (i in 1:NROW(df)){
        temp_mat[i,] <- boot.ci(boot_obj,
                                type = "bca",
                                conf = 1 - alpha,
                                index = i)$bca[4:5]
    }

    lwr <- temp_mat[,1]
    upr <- temp_mat[,2]

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- out
    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}
