## add_probs method for lm objects
## TODO : Use t or Z distribution option
## TODO : change "Y" to the name of the response

add_probs.lm <- function(tb, fit, quant, probName = NULL,
                         direction = "<"){
    if (is.null(probName) && direction == "<")
        probName <- paste("Pr(Y < ", quant, ")", sep="")
    if (is.null(probName) && direction == ">")
        probName <- paste("Pr(Y > ", quant, ")", sep="")
    if ((probName %in% colnames(tb))) {
        warning ("These probabilities may have already been appended to your dataframe")
        return(tb)
    }
    out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
    fitted <- out$fit[,1]
    residual_df <- out$df
    se_fitted <- out$se.fit
    resid_var <- out$residual.scale^2
    se_pred <- sqrt(resid_var + se_fitted^2)
    t_quantile <- (quant - fitted) / se_pred
    if (direction == "<")
        t_prob <- pt(q = t_quantile, df = residual_df)
    if (direction == ">")
        t_prob <- 1 - pt(q = t_quantile, df = residual_df)
    if (is.null(tb[["pred"]]))
        tb[["pred"]] <- fitted
    if (is.null(tb[[probName]]))
        tb[[probName]] <- t_prob
    return(tb)
}
