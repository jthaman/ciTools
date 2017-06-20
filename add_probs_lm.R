## add_probs method for lm objects
add_probs.lm <- function(tb, fit, quantile = 0.20, probName = NULL){
    if (quant <= 0 || quant >= 1)
        stop ("quant should be in (0,1)")
    if (is.null(quantileName))
        quantileNameCheck <- paste(quant, "-Quantile", sep="")
    if (!(quantileNameCheck %in% colnames(tb)))
        quantileName <- quantileNameCheck
    out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
    fitted <- out$fit[,1]
    residual_df <- out$df
    se_fitted <- out$se.fit
    t_quantile <- qt(p = quant, df = residual_df)
    out_quantiles <- fitted + se_fitted * t_quantile
    if (is.null(tb[["pred"]])) tb[["pred"]] <- fitted
    tb[[quantileName]] <- out_quantiles
    tb
    
}


