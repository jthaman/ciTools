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

#' Confidence Intervals for Generalized Linear Mixed Model Predictions
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{glmerMod}. Confidence intervals are approximate and
#' determined via simulation.
#'
#' There are two methods one can use to calculate confidence intervals
#' for GLMM fits: bootstrap or parametric. The default and recommended
#' method is bootstrap. The bootstrap method can handle many types of
#' models and we find it to be generally reliable and robust as it is
#' built on the \code{bootMer} function from \code{lme4}. An
#' experimental parametric method is included that mimics the
#' functionality of \code{add_ci.lmer}'s default method. We caution
#' against using this method because presently it only works for GLMMs
#' that have a single random intercept term.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{glmerMod}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds automatically will be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param type A string. If \code{type == "boot"} then bootstrap
#'     intervals are formed. If \code{type == "parametric"} then
#'     parametric intervals are formed.
#' @param yhatName \code{NULL} or a string. Name of the predictions
#'     vector. If \code{NULL}, the predictions will be named
#'     \code{pred}.
#' @param response A logical. The default is \code{TRUE}. If
#'     \code{TRUE}, the confidence intervals will be determined for
#'     the expected response; if \code{FALSE}, confidence intervals
#'     will be made on the scale of the linear predictor.
#' @param includeRanef A logical. Default is \code{TRUE}. Set whether
#'     the predictions and intervals should be made conditional on the
#'     random effects. If \code{FALSE}, random effects will not be
#'     included.
#' @param nSims A positive integer.  Controls the number of bootstrap
#'     replicates if \code{type = "boot"}.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @seealso \code{\link{add_pi.glmerMod}} for prediction intervals
#'     of \code{glmerMod} objects, \code{\link{add_probs.glmerMod}} for
#'     conditional probabilities of \code{glmerMod} objects, and
#'     \code{\link{add_quantile.glmerMod}} for response quantiles of
#'     \code{glmerMod} objects.
#'
#' @references
#' http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
#'
#' @examples
#' ## random intercept example
#' tb <- data.frame(y=rpois(1000,lambda=3),x=runif(1000),
#'                  f=factor(sample(1:10,size=1000,replace=TRUE)))
#' fit <- lme4::glmer(y~x+(1|f),data=tb,family=poisson)
#'
#' add_ci(tb, fit, includeRanef = TRUE, names = c("LCB", "UCB"), type = "parametric")
#' \dontrun{add_ci(tb, fit, includeRanef = TRUE,
#'                 names = c("LCBB", "UCBB"), type = "boot")}
#'
#' @export

add_ci.glmerMod <- function(tb, fit,
                            alpha = 0.05, names = NULL, yhatName = "pred",
                            response = TRUE,
                            type = "boot", includeRanef = TRUE,
                            nSims = 250, ...){

    if (!is.null(fit@optinfo$conv$lme4$code))
        warning ("Coverage probabilities may be inaccurate if the model failed to converge")

    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe. Overwriting.")
    }

    if (type == "boot") {
        bootstrap_ci_glmermod(tb, fit, alpha, names, includeRanef, nSims, yhatName, response)
    } else if (type == "parametric") {
        parametric_ci_glmermod(tb, fit, alpha, names, includeRanef, yhatName, response)
    } else {
        stop("Incorrect type specified!")
    }
}

parametric_ci_glmermod <- function(tb, fit, alpha, names, includeRanef, yhatName, response){

    if (length(fit@cnms[[1]]) != 1)
        stop("parametric confidence intervals are currently only implemented for random intercept models.")

    seFixed <- get_prediction_se_mermod(tb, fit)
    ranef_name <- names(fit@cnms)[1] ## just one random effect for now
    seRandom <- arm::se.ranef(fit)[[1]][,1]
    seRandom_vec <- rep(NA, length(tb[[ranef_name]]))

    seRandom_df <- tibble::tibble(
                               group = names(seRandom),
                               seRandom = seRandom
    )

    names(seRandom_df)[names(seRandom_df) == 'group'] <- ranef_name
    seRandom_vec <- dplyr::left_join(tb, seRandom_df, by = ranef_name)[["seRandom"]]

    rdf <- get_resid_df_mermod(fit)

    if (fit@resp$family$family %in% c("binomial", "poisson")){
        crit_val <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1)
    } else {
        crit_val <- qt(p = 1 - alpha/2, df = rdf)
    }

    inverselink <- fit@resp$family$linkinv

    if(includeRanef) {
        re.form <- NULL
        seGlobal <- sqrt(seFixed^2 + seRandom_vec^2)
    } else {
        re.form <-  NA
        seGlobal <- seFixed
    }

    out <- predict(fit, tb, re.form = re.form)

    pred <- out
    upr <- out + crit_val * seGlobal
    lwr <- out - crit_val * seGlobal

    if (response == TRUE){
        pred <- inverselink(pred)
        upr <- inverselink(upr)
        lwr <- inverselink(lwr)
    }

    if(fit@resp$family$link %in% c("inverse", "1/mu^2")){
        upr1 <- lwr
        lwr <- upr
        upr <- upr1
    }

    tb[[yhatName]] <- pred
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)
}

ciTools_data <- new.env(parent = emptyenv())

bootstrap_ci_glmermod <- function(tb, fit, alpha, names, includeRanef, nSims, yhatName, response) {

    ciTools_data$tb_temp <- tb

    if (includeRanef) {
        rform = NULL
        my_pred <- my_pred_full_glmer
    } else {
        rform = NA
        my_pred <- my_pred_fixed_glmer
    }

    if (response){
        lvl <- "response"
    }
    else{
        lvl <- "link"
    }

    boot_obj <- lme4::bootMer(fit, my_pred, nsim=nSims, type="parametric", re.form = rform, lvl = lvl)
    ci_out <- boot_quants(boot_obj, alpha)

    tb[[yhatName]] <- predict(fit, tb, re.form = rform, type = lvl)
    tb[[names[1]]] <- ci_out$lwr
    tb[[names[2]]] <- ci_out$upr
    tibble::as_data_frame(tb)
}
