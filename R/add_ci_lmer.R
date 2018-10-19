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

#' Confidence Intervals for Linear Mixed Model Predictions
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{lmerMod}. It is recommended that one use parametric
#' confidence intervals when modeling with a random intercept linear
#' mixed model (i.e. a fit with a formula such as \code{lmer(y ~ x +
#' (1|group))}). Otherwise, confidence intervals may be bootstrapped
#' via \code{lme4::bootMer}.
#'
#' Bootstrapped intervals are slower to compute, but they are the
#' recommended method when working with any linear mixed models more
#' complicated than the random intercept model.
#'
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lmerMod}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, confidence bounds automatically will be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param type A string. Must be \code{"parametric"} or \code{"boot"},
#'     If \code{type = "boot"}, then \code{add_ci} calls
#'     \code{lme4::bootMer} to calculate the confidence
#'     intervals. This method may be time consuming, but is applicable
#'     with random slope and random intercept models. The parametric
#'     method is fast, but currently only works well for random
#'     intercept models.
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
#' @seealso \code{\link{add_pi.lmerMod}} for prediction intervals
#'     of \code{lmerMod} objects, \code{\link{add_probs.lmerMod}} for
#'     conditional probabilities of \code{lmerMod} objects, and
#'     \code{\link{add_quantile.lmerMod}} for response quantiles of
#'     \code{lmerMod} objects.
#'
#' @examples
#' \dontrun{
#' dat <- lme4::sleepstudy
#' # Fit a linear mixed model (random intercept model)
#' fit <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#' # Get the fitted values for each observation in dat, and
#' # append CIs for those fitted values to dat
#' add_ci(dat, fit, alpha = 0.5)
#' # Try the parametric bootstrap method, and make prediction at the population level
#' add_ci(dat, fit, alpha = 0.5, type = "boot", includeRanef = FALSE, nSims = 100)
#' }
#' @export

add_ci.lmerMod <- function(tb, fit,
                           alpha = 0.05, names = NULL, yhatName = "pred",
                           type = "boot", includeRanef = TRUE,
                           nSims = 500, ...){

  if(!is.null(attr(fit@pp$X, "msgRankdrop")))
    warning("Model matrix is rank deficient!")

  if (!is.null(fit@optinfo$conv$lme4$code))
    warning ("Coverage probabilities may be inaccurate if the model failed to converge")

  if (is.null(names)){
    names[1] <- paste("LCB", alpha/2, sep = "")
    names[2] <- paste("UCB", 1 - alpha/2, sep = "")
  }
  if ((names[1] %in% colnames(tb))) {
    warning ("These CIs may have already been appended to your dataframe. Overwriting.")
  }

  if (type == "parametric")
    parametric_ci_lmermod(tb, fit, alpha, names, includeRanef, yhatName)
  else if (type == "boot")
    bootstrap_ci_lmermod(tb, fit, alpha, names, includeRanef, nSims, yhatName)
  else
    stop("Incorrect type specified!")
}


parametric_ci_lmermod <- function(tb, fit, alpha, names, includeRanef, yhatName){

  if (length(fit@cnms[[1]]) != 1)
    stop("parametric confidence intervals are currently only implemented for random intercept models.")

  seFixed <- get_prediction_se_mermod(tb, fit)
  seRandom <- arm::se.ranef(fit)[[1]][1,]

  rdf <- get_resid_df_mermod(fit)

  if(includeRanef) {
    re.form <- NULL
    seGlobal <- sqrt(seFixed^2 + seRandom^2)
  } else {
    re.form <-  NA
    seGlobal <- seFixed
  }

  out <- predict(fit, tb, re.form = re.form)
  if(is.null(tb[[yhatName]]))
    tb[[yhatName]] <- out
  tb[[names[1]]] <- out + qt(alpha/2, df = rdf) * seGlobal
  tb[[names[2]]] <- out + qt(1 - alpha/2, df = rdf) * seGlobal
  tibble::as_data_frame(tb)

}


ciTools_data <- new.env(parent = emptyenv())

bootstrap_ci_lmermod <- function(tb, fit, alpha, names, includeRanef, nSims, yhatName) {

  ciTools_data$tb_temp <- tb

  if (includeRanef) {
    rform = NULL
    my_pred <- my_pred_full
  } else {
    rform = NA
    my_pred <- my_pred_fixed
  }

  boot_obj <- lme4::bootMer(fit, my_pred, nsim=nSims, type="parametric", re.form = rform)

  ci_out <- boot_quants(boot_obj, alpha)

  if(is.null(tb[[yhatName]]))
    tb[[yhatName]] <- ci_out$fit
  tb[[names[1]]] <- ci_out$lwr
  tb[[names[2]]] <- ci_out$upr
  as_data_frame(tb)

}
