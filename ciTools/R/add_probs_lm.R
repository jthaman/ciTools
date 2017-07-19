#' Event Probabilities for Linear Models
#'
#' This is the method \code{add_probs} uses if the model fit is
#' linear.  Probabilities are calculated using the t-distribution.
#' 
#'
#' @param tb A tibble or Data Frame on which to append probabilities
#' @param fit An object of class lm. Predictions are made with this
#'     object.
#' @param q A double. A quantile of the response distribution.
#' @param name NULL or character vector of length one. If \code{NULL},
#'     probabilities will automatically be named by
#'     \code{add_probs()}, otherwise, the probabilities will be named
#'     \code{name} in the returned tibble
#' @param comparison A character vector of length one. If
#'     \code{comparison = "<"}, then Pr(Y|x < q) is calculated for
#'     each observation in \code{tb}. Must be "<" or ">" for linear,
#'     log-linear and linear mixed models. If \code{fit} is a glm,
#'     then \code{comparison} may also be "<=", ">=", or "=".
#' @param log_response A logical. If TRUE, then the quantile
#'     \code{t} will be transformed to the log scale.
#' @return A tibble, \code{tb}, with predicted values and
#'     probabilities attached.
#' 
#' @examples
#' # linear regression
#' fit1 <- lm(dist ~ speed, data = cars)
#' # Calculate Pr(dist < 40) for each observation in cars
#' add_probs(cars, fit1, q = 40)
#' add_probs(cars, fit1, q = 50, comparison = ">")
#' 
#' @export

add_probs.lm <- function(tb, fit, q, name = NULL,
                         comparison = "<", log_response = FALSE){

    if (is.null(name) && comparison == "<")
        name <- paste("prob_less_than", q, sep="")
    if (is.null(name) && comparison == ">")
        name <- paste("prob_greater_than", q, sep="")
    if ((name %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe. Overwriting.")
    }

    if (log_response)
        q <- log(q)

    out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
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
    if (is.null(tb[["pred"]]))
        tb[["pred"]] <- fitted
    if (is.null(tb[[name]]))
        tb[[name]] <- t_prob
    tibble::as_data_frame(tb)
}
