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

#' Quantiles for the Response of a Linear Mixed Model
#'
#' This function is one of the methods for \code{add_quantile} and
#' is called automatically when \code{add_quantile} is applied to an
#' object of class \code{lmerMod}.
#'
#' \code{add_qauntile.lmerMod} may use one of two different methods
#' for determining quantiles: a parametric method or a parametric
#' bootstrap method (via \code{lme4::simulate}). The parametric method
#' is the default. Only use the parametric method (\code{type =
#' "parametric"}) if \code{fit} is a random intercept model,
#' e.g. \code{fit = lmer(y ~ x + (1|group))}. If your model of
#' interest is random slope and random intercept, use the parametric
#' bootstrap method (\code{type = "boot"}).
#' 
#' @param tb A tibble or data frame of new data.
#' @param fit An object of class \code{lm}. Predictions are made with
#'     this object.
#' @param p A real number between 0 and 1.  Determines the probability
#'     level of the quantiles.
#' @param name \code{NULL} or a string. If \code{NULL}, quantiles
#'     automatically will be named by \code{add_quantile}, otherwise,
#'     they will be named \code{name}.
#' @param yhatName A string. Determines the name of the column of
#'     predictions.
#' @param includeRanef The random effects be included or not? If
#'     \code{TRUE}, quantiles will be calculated at the
#'     "group level". Otherwise, quantiles will be calculated at the
#'     "population level", where random effects are set to \eqn{0}.
#' @param type A string. Options are \code{"parametric"} or
#'     \code{"boot"}.
#' @param nSims A positive integer. Set the number of bootstrap
#'     simulations to perform. Only applied when \code{type = "boot"}.
#' @param log_response A logical. Set to \code{TRUE} if the model is a
#'     log-linear mixed model: \eqn{\log(Y) = X\beta + Z\gamma +
#'     \epsilon}.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values and level-p
#'     quantiles attached.
#' 
#' @seealso \code{\link{add_ci.lmerMod}} for confidence intervals
#'     for \code{lmerMod} objects, \code{\link{add_pi.lmerMod}} for
#'     prediction intervals of \code{lmerMod} objects, and
#'     \code{\link{add_probs.lmerMod}} for response probabilities of
#'     \code{lmerMod} objects.
#'
#' @examples
#' dat <- lme4::sleepstudy
#'
#' # Fit a random intercept model
#' fit <- lme4::lmer(Reaction ~ Days + (1|Subject), data = lme4::sleepstudy)
#'
#' # Using the parametric method: given the model fit, what value
#' # of reaction time do we expect half of new reaction times to fall
#' # under?
#' add_quantile(dat, fit, p = 0.5)
#'
#' # Using the parametric method:
#' # as above, but we ignore the random effects.
#' add_quantile(dat, fit, p = 0.5, includeRanef = FALSE)
#'
#' @export

add_quantile.lmerMod <- function(tb, fit, 
                                 p, name = NULL,
                                 yhatName = "pred", includeRanef = TRUE,
                                 type = "boot",
                                 nSims = 200, log_response = FALSE, ...) {

    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)){
        name <- paste("quantile", p, sep = "")
    }
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }
    if (type == "parametric")
        parametric_quantile_lmermod(tb, fit, p, name, includeRanef, log_response, yhatName)
    else if (type == "boot")
        boot_quantile_lmermod(tb, fit, p, name, includeRanef, nSims, log_response, yhatName)
    else
        stop ("Incorrect type specified")
}

parametric_quantile_lmermod <- function(tb, fit, p, name, includeRanef, log_response, yhatName){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    out <- predict(fit, tb, re.form = re.form)
    quant <- out + qt(p ,df = rdf) * seGlobal
    if (log_response)
        out <- exp(out)
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quant
    if(log_response)
        tb[[name]]<- exp(quant)
    tibble::as_data_frame(tb)
}


boot_quantile_lmermod <- function(tb, fit, p, name, includeRanef, nSims, log_response, yhatName) {

    if (includeRanef) 
        reform = NULL
    else 
        reform = NA

    gg <- simulate(fit, newdata = tb, re.form = reform, nsim = nSims)
    gg <- as.matrix(gg)
    quant <- apply(gg, 1, FUN = quantile, probs = p)

    out <- predict(fit, tb, re.form = reform)
    if (log_response)
        out <- exp(out)
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- out
    tb[[name]] <- quant
    if (log_response)
        tb[[name]]<- exp(quant)
    tibble::as_data_frame(tb)
}
