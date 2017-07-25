#' Confidence Intervals for Expected Response of Generalized Linear Models
#'
#' This function is one of the methods for \code{add_ci}, and is
#' called automatically when \code{add_ci} is used on a \code{fit} of
#' class \code{glm}. Confidence Intervals are determined by making an
#' interval on the scale of the linear predictor, then applying the
#' inverse link function from the model fit.
#'
#' @param tb A tibble or Data Frame.
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names NULL or character vector of length two. If
#'     \code{NULL}, confidence bounds will automatically be named by
#'     \code{add_ci}, otherwise, the lower confidence bound will be
#'     named \code{names[1]} and the upper confidence bound will be
#'     named \code{names[2]}.
#' @param type A string, either \code{"parametric"} or
#'     \code{"bootstrap"}.  but at the moment, the bootstrap method is
#'     not implemented.
#' @param response A logical. If \code{TRUE}, the confidence intervals
#'     will be determined for the expected response, if \code{FALSE},
#'     confidence intervals will be made on the scale of the linear
#'     predictor.
#' @param ... Additional arguments.
#' @return A tibble, \code{tb}, with predicted values, upper and lower
#'     confidence bounds attached.
#'
#' 
#' @export

add_ci.glm <- function(tb, fit, alpha = 0.05, names = NULL,
                      response = TRUE, type = "parametric", ...){

    if (grepl("numerically 0 or 1", list(warnings())))
        warning ("If there is perfect separation in your logistic regression, you shouldn't trust these confidence intervals")
    if (is.null(names)){
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe. Overwriting.")
    }
    if (type == "bootstrap")
        stop ("not yet implemented")
    else if (type == "parametric")
        parametric_ci_glm(tb, fit, alpha, names, response)
    else
        if(!(type %in% c("bootstrap", "parametric"))) stop("Incorrect interval type specified!")

}

parametric_ci_glm <- function(tb, fit, alpha, names, response){
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
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- pred
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr
    tibble::as_data_frame(tb)

}


## parametric_ci_glm <- function(tb, fit, alpha, names, type = "response"){

##     inverselink <- fit$family$linkinv
##     if(type == "response"){
##         out <- glm_ci_response(tb, fit, alpha, inverselink) %>%
##             plyr::rename(c("lwr" = names[1], "upr" = names[2])) %>%
##             select(-.se)
##     }

##     if(type == "link"){
##         out <- glm_ci_linear_predictor(tb, fit, alpha) %>%
##             plyr::rename(c("lwr" = names[1], "upr" = names[2])) %>%
##             select(-.se)
##     }
##     if(!(type %in% c("response", "link"))) top("Incorrect interval type specified!")
##     out
## }

## glm_ci_response <- function(tb, fit, alpha, ilink){

##     if(is.null(tb[["pred"]]))
##         tb <- modelr::add_predictions(tb, fit)
##     tb %>%
##         add_standard_error(fit) %>%
##         mutate(
##             lwr = ilink(pred + qt(alpha/2, df = fit$df.residual) * .se),
##             upr = ilink(pred + qt(1 - alpha/2, df = fit$df.residual) * .se)
##         )

## }

## glm_ci_linear_predictor <- function(tb, fit, alpha){

##     if(is.null(tb[["pred"]]))
##         tb <- modelr::add_predictions(tb, fit)
##     tb %>%
##         add_standard_error(fit) %>%
##         mutate(
##             lwr = pred + qt(alpha/2, df = fit$df.residual) * .se,
##             upr = pred + qt(1 - alpha/2, df = fit$df.residual) * .se
##         )

## }
