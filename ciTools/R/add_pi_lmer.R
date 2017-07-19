#' Prediction Intervals for the Response of Linear Mixed Models
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{lmerMod}. It is recommended that one use parametric
#' prediction intervals when modeling with a random intercept
#' LMM. Otherwise prediction intervals may be simulated using one of
#' two methods. \code{"sim"} indicates that \code{predictInterval}
#' from \code{merTools} should be used to simulate responses to form
#' prediction intervals, and \code{"sim_lme4"} indicates that
#' \code{simulate.merMod} should be used to simulate predictions.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class \code{lmerMod}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, prediction bounds will automatically be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param type A string, either \code{"parametric"}, \code{"sim"},
#'     \code{"sim_lme4"}.
#' @param includeRanef A logical. Set whether the predictions and
#'     intervals should be made conditional on the random effects. If
#'     \code{FALSE}, random effects will not be included.
#' @param nSims A positive integer. If \code{type = "sim"} or
#'     \code{type = "sim_lme4"}, \code{nSims} will determine the
#'     number of simulated draws to make.
#' @param log_response A logical, indicating if the response is on log
#'     scale in the model. If \code{TRUE}, prediction intervals will
#'     be returned on the response scale.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @examples
#' \dontrun{
#' data(sleepstudy) ## included in lme4
#' tb <- sleepstudy
#' fit1 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
#' add_pi.lmerMod(tb, fit1)
#' add_pi.lmerMod(tb, fit1, type = "sim")
#' add_pi.lmerMod(tb, fit1, type = "sim_lme4")
#' }
#' 
#' @export

add_pi.lmerMod <- function(tb, fit, 
                          alpha = 0.05, type = "parametric", includeRanef = TRUE,
                          names = NULL, log_response = FALSE, nSims = 200, ...) {
    if (is.null(names)){
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    }

    if(type == "sim") 
        sim_pi_mermod(tb, fit, alpha, names, includeRanef, nSims, log_response, ...)
    else if(type == "parametric")
        parametric_pi_mermod(tb, fit, alpha, names, includeRanef, log_response, ...)
    else if(type == "sim_lme4")
        lme4_pi_mermod(tb, fit, alpha, names, includeRanef, nSims, log_response, ...)
    else
        stop("Incorrect type specified!")

 }

sim_pi_mermod <- function(tb, fit, alpha, names, includeRanef, nSims, log_response) {

    if (includeRanef) {
        which = "full"
        reform = NULL
    } else {
        which = "fixed"
        reform = NA
    }

    pi_out <- suppressWarnings(predictInterval(fit, tb, which = which, level = 1 - alpha,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = TRUE))

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = reform)
    tb[[names[1]]] <- pi_out$lwr
    tb[[names[2]]] <- pi_out$upr
    if (log_response){
        tb[["pred"]] <- exp(predict(fit, tb, re.form = reform))
        tb[[names[1]]] <- exp(tb[[names[1]]])
        tb[[names[2]]] <- exp(tb[[names[2]]])
    }

    tibble::as_data_frame(tb)
    
}

parametric_pi_mermod <- function(tb, fit, alpha, names, includeRanef, log_response){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    out <- predict(fit, tb, re.form = re.form)
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[names[1]]] <- out + qt(alpha/2,df = rdf) * seGlobal
    tb[[names[2]]] <- out + qt(1 - alpha/2, df = rdf) * seGlobal
    if (log_response){
        tb[["pred"]] <- exp(out)
        tb[[names[1]]] <- exp(tb[[names[1]]])
        tb[[names[2]]] <- exp(tb[[names[2]]])
    }
    tibble::as_data_frame(tb)
}



lme4_pi_mermod <- function(tb, fit, alpha, names, includeRanef, nSims, log_response) {

    if (includeRanef) 
        reform = NULL
    else 
        reform = NA

    gg <- simulate(fit, re.form = reform, nsim = nSims)
    gg <- as.matrix(gg)
    lwr <- apply(gg, 1, FUN = quantile, probs = alpha/2)
    upr <- apply(gg, 1, FUN = quantile, probs = 1 - alpha / 2)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = reform)
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    if (log_response){
        tb[["pred"]] <- exp(predict(fit, tb, re.form = reform))
        tb[[names[1]]] <- exp(tb[[names[1]]])
        tb[[names[2]]] <- exp(tb[[names[2]]])
    }
    tibble::as_data_frame(tb)
}

## parametric_pi_mermod <- function(tb, fit, alpha, names, includeRanef){
##     if (includeRanef == TRUE)
##         reform = NULL
##     else
##         reform = NA

##     X <- model.matrix(reformulate(attributes(terms(fit))$term.labels), tb)
##     vcovBetaHat <- vcov(fit)
    
##     seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
##         diag() %>%
##         sqrt()
    
##     seRandom <- arm::se.ranef(fit)[[1]][1,]
##     rdf <- nrow(model.matrix(fit)) - length(fixef(fit)) -
##         (length(attributes(summary(fit)$varcor)$names) + 1)
##     se_residual <- sigma(fit)
##     if(includeRanef)
##         seGlobal <- sqrt(seFixed^2 + seRandom^2 + se_residual^2)
##     else
##         seGlobal <- sqrt(seFixed^2 + se_residual^2)
##     if(is.null(tb[["pred"]]))
##         tb[["pred"]] <- predict(fit, tb, re.form = reform) 
##     tb[[names[1]]] <- tb[["pred"]] + qt(alpha/2, df = rdf) * seGlobal
##     tb[[names[2]]] <- tb[["pred"]] + qt(1 - alpha/2, df = rdf) * seGlobal
##     tb
    
## }



## method that uses simulate from lme4

