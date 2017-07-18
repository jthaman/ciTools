#' Confidence Intervals for the Response of Linear Mixed Models
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{lmerMod}. It is recommended that one use parametric
#' confidence intervals when modeling with a random intercept
#' LMM. Otherwise confidence intervals may be simulated (type =
#' \code{"sim"}) via \code{predictInterval} from \code{merTools}. A
#' bootstrap method may be included in the future.
#' 
#' @param tb A tibble or Data Frame.
#' @param fit An object of class \code{lmerMod}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param type A string, either \code{"parametric"} or
#'     \code{"bootstrap"}. But at the moment, the bootstrap method is
#'     not implemented.
#' @param includeRanef A logical. Set whether the predictions and
#'     intervals should be made conditional on the random effects. If
#'     \code{FALSE}, random effects will not be included.
#' @param nSims A positive integer. If \code{type = "sim"},
#'     \code{nSims} will determine the number of simulated draws to
#'     make.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' @examples
#' \dontrun{
#' data(sleepstudy) ## included in lme4
#' tb <- sleepstudy
#' fit1 <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
#' add_ci(tb, fit1)
#' add_ci(tb, fit1, type = "sim")
#' }
#' 
#' @export

add_ci.lmerMod <- function(tb, fit, 
                           alpha = 0.05, type = "parametric", includeRanef = TRUE,
                           names = NULL, nSims = 200, ...){

    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe")
        return(tb)
    }

    if (type == "parametric")
        parametric_ci_mermod(tb, fit, alpha, names, includeRanef)
    else if (type == "sim")
        sim_ci_mermod(tb, fit, alpha, names, includeRanef, nSims)
    else if (type == "boot")
        bootstrap_ci_mermod(tb, fit, alpha, names, includeRanef, nSims)
    else
        stop("Incorrect type specified!")
}


parametric_ci_mermod <- function(tb, fit, alpha, names, includeRanef){
    
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
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out
    tb[[names[1]]] <- out + qt(alpha/2, df = rdf) * seGlobal
    tb[[names[2]]] <- out + qt(1 - alpha/2, df = rdf) * seGlobal
    tibble::as_data_frame(tb)
    
}


## parametric_ci_mermod <- function(tb, fit, alpha, names, includeRanef){
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
##     if(includeRanef)
##         seGlobal <- sqrt(seFixed^2 + seRandom^2)
##     else
##         seGlobal <- seFixed
##     if(is.null(tb[["pred"]]))
##         tb[["pred"]] <- predict(fit, tb, re.form = reform) 
##     tb[[names[1]]] <- tb[["pred"]] + qt(alpha/2, df = rdf) * seGlobal
##     tb[[names[2]]] <- tb[["pred"]] + qt(1 - alpha/2, df = rdf) * seGlobal
##     tb
## }

## simulation method using predictInterval
sim_ci_mermod <- function(tb, fit, alpha, names, includeRanef, nSims = 200) {

    if (includeRanef) {
        which = "full"
        reform = NULL
    } else {
        which = "fixed"
        reform = NA
    }

    ci_out <- suppressWarnings(predictInterval(fit, tb, which = which, level = 1 - alpha,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = FALSE))

    if(is.null(tb[["pred"]])) 
        tb[["pred"]] <- predict(fit, tb, re.form = reform) 
    tb[[names[1]]] <- ci_out$lwr
    tb[[names[2]]] <- ci_out$upr
    tibble::as_data_frame(tb)
    
}

bootstrap_ci_mermod <- function(tb, fit, alpha, names, includeRanef, nSims) {
    
    if (includeRanef) { 
        rform = NULL
        my_pred <- my_pred_full
    } else {
        rform = NA
        my_pred <- my_pred_fixed
    }
        
    boot_obj <- lme4::bootMer(fit, my_pred, nsim=nSims, type="parametric", re.form = rform)

    ci_out <- boot_quants(boot_obj, alpha) 

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- ci_out$fit
    tb[[names[1]]] <- ci_out$lwr
    tb[[names[2]]] <- ci_out$upr
    tb
    
}

## ##Bootstrap CIs for merMod objects
## ##Note that the bootstrapping is done on the random effects only. 

## oldbootstrap_ci_mermod <- function(tb, fit, alpha, names, nBS = 999){
    
##     X <- model.matrix(fit)
##     vcovBetaHat <- vcov(fit)
    
##     seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
##         diag() %>%
##         sqrt()
    
##     wholePlotBS <- sample_ranefs(fit, nBS, alpha) 
##     if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
##     tb[[names[1]]] <- tb[["pred"]] + qnorm(alpha/2) * seFixed + wholePlotBS[1]
##     tb[[names[2]]] <- tb[["pred"]] + qnorm(1 - alpha/2) * seFixed + wholePlotBS[2]
##     tb
    
## }

## sample_ranefs <- function(fit, nBS, alpha){
##     if(length(ranef(fit)[[1]][[1]]) < 5) warning(
##                                              "Fewer than 5 whole plots; consider using parametric methods")
##     ranef(fit)[[1]][[1]] %>% 
##         base::sample(size = nBS, replace = T) %>%
##         quantile(probs = c(alpha/2, 1-alpha/2))
## }
