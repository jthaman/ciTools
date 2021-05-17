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

#' Confidence Intervals for Generalized Linear Model Predictions
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{glm}. The default method calculates confidence
#' intervals by making an interval on the scale of the linear
#' predictor, then applying the inverse link function from the model
#' fit to transform the linear level confidence intervals to the
#' response level. Alternatively, confidence intervals may be
#' calculated through a nonparametric bootstrap method.
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds automatically will be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param yhatName A character vector of length one. Name of the
#'     vector of predictions made for each observation in df
#' @param type A character vector of length one. Must be \code{type =
#'     "parametric"} or \code{type = "boot"}. \code{type} determines
#'     the method used to compute the confidence intervals.
#' @param response A logical. The default is \code{TRUE}. If
#'     \code{TRUE}, the confidence intervals will be determined for
#'     the expected response; if \code{FALSE}, confidence intervals
#'     will be made on the scale of the linear predictor.
#' @param nSims An integer. Number of simulations to perform if the
#'     bootstrap method is used.
#' @param ... Additional arguments passed to \code{boot::boot()}.
#'
#' @return A dataframe, \code{df}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @seealso \code{\link{add_pi.glm}} for prediction intervals for
#'     \code{glm} objects, \code{\link{add_probs.glm}} for conditional
#'     probabilities of \code{glm} objects, and
#'     \code{\link{add_quantile.glm}} for response quantiles of
#'     \code{glm} objects.
#'
#' @examples
#' # Poisson regression
#' fit <- glm(dist ~ speed, data = cars, family = "poisson")
#' add_ci(cars, fit)
#' # Try a different confidence level
#' add_ci(cars, fit, alpha = 0.5)
#' # Add custom names to the confidence bounds (may be useful for plotting)
#' add_ci(cars, fit, alpha = 0.5, names = c("lwr", "upr"))
#'
#' # Logistic regression
#' fit2 <- glm(I(dist > 30) ~ speed, data = cars, family = "binomial")
#' dat <- cbind(cars, I(cars$dist > 30))
#' # Form 95% confidence intervals for the fit:
#' add_ci(dat, fit2)
#' # Form 50% confidence intervals for the fit:
#' add_ci(dat, fit2, alpha = 0.5)
#' # Make confidence intervals on the scale of the linear predictor
#' add_ci(dat, fit2, alpha = 0.5, response = FALSE)
#' # Add custom names to the confidence bounds
#' add_ci(dat, fit2, alpha = 0.5, names = c("lwr", "upr"))
#'
#'
#' @export

add_ci.glm <- function(df, fit, alpha = 0.05, names = NULL, yhatName = "pred",
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
        boot_ci_glm(
          df=df,
          fit=fit,
          alpha=alpha,
          names=names,
          yhatName=yhatName,
          response=response,
          nSims=nSims,
          ...
        )
    else if (type == "parametric")
        parametric_ci_glm(df, fit, alpha, names, yhatName, response)
    else
        stop("Incorrect interval type specified!")

}

parametric_ci_glm <- function(df, fit, alpha, names, yhatName, response){
    out <- predict(fit, df, se.fit = TRUE, type = "link")

    if (fit$family$family %in% c("binomial", "poisson"))
        crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
    else
        crit_val <- qt(p = 1 - alpha/2, df = fit$df.residual)

    inverselink <- fit$family$linkinv

    if (response){
        pred <- inverselink(out$fit)
        upr <- inverselink(out$fit + crit_val * out$se.fit)
        lwr <- inverselink(out$fit - crit_val * out$se.fit)
        if(fit$family$link %in% c("inverse", "1/mu^2")){
            ## need to do something like this for any decreasing link
            ## function.
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

boot_fit <- function(data, df, fit, lvl, indices){
    data_temp <- data[indices,]
    form <- fit$formula
    fam <- fit$family
    ## temp_fit <- glm(form, data = data_temp, family = fam)
    temp_fit <- update(fit, data = data_temp)
    # Find factor levels that are in the prediction data that are not in the
    # updated fit data
    factor_columns <-
      names(temp_fit$model[, sapply(X=temp_fit$model, FUN=is.factor), drop=FALSE])
    missing_factors <-
      lapply(
        X=setNames(nm=factor_columns),
        FUN=function(current_nm) {
          base::setdiff(
            unique(df[[current_nm]]),
            unique(data_temp[[current_nm]])
          )
        }
      )
    # Set the parameter values to NA when they are not in the model
    df_new <-
      lapply(
        X=setNames(nm=names(df)),
        FUN=function(current_nm) {
          ret <- df[[current_nm]]
          mask_na <- ret %in% missing_factors[[current_nm]]
          if (any(mask_na)) {
            warning(
              "Factor levels in column ", current_nm, " for prediction are missing in the model.  ",
              "Consider setting `strata` to ensure that all factors are represented."
            )
            ret[mask_na] <- NA
          }
          ret
        }
      )
    predict(temp_fit, newdata = as.data.frame(df_new), type = lvl)
}

boot_ci_glm <- function(df, fit, alpha, names, yhatName, response, nSims, ...){
    if (response){
        lvl <- "response"
    }
    else{
        lvl <- "link"
    }

    out <- predict(fit, df, type = lvl)

    boot_obj <- boot(data = fit$model,
                     statistic = boot_fit,
                     R = nSims,
                     fit = fit,
                     lvl = lvl,
                     df = df,
                     ...)

    temp_mat <- matrix(0, ncol = 2, nrow = NROW(df))

    for (i in 1:NROW(df)){
        temp_mat[i,] <- boot.ci(boot_obj,
                                type = "bca",
                                conf = 1 - alpha,
                                index = i)$bca[4:5]
    }

    lwr <- temp_mat[,1]
    upr <- temp_mat[,2]

    ## raw_boot <- boot_obj$t
    ## lwr <- apply(raw_boot, 2, FUN = quantile, probs = alpha / 2)
    ## upr <- apply(raw_boot, 2, FUN = quantile, probs = 1 - alpha / 2)

    if(is.null(df[[yhatName]]))
        df[[yhatName]] <- out
    df[[names[1]]] <- lwr
    df[[names[2]]] <- upr
    data.frame(df)
}
