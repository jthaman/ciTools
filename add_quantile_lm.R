## add_quantile method for lm objects

add_quantile.lm <- function(tb, fit, prob = 0.20, quantileName = NULL){
    if (prob == 0.5)
        warning ("The 0.5 quantile is equal to the fitted values")
    if (prob <= 0 || prob >= 1)
        stop ("prob should be in (0,1)")
    if (is.null(quantileName))
        quantileNameCheck <- paste(prob, "-quantile", sep="")
    if (!(quantileNameCheck %in% colnames(tb)))
        quantileName <- quantileNameCheck
    out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
    fitted <- out$fit[,1]
    residual_df <- out$df
    se_fitted <- out$se.fit
    t_quantile <- qt(p = prob, df = residual_df)
    out_quantiles <- fitted + se_fitted * t_quantile
    if (is.null(tb[["pred"]])) tb[["pred"]] <- fitted
    tb[[quantileName]] <- out_quantiles
    return(tb)
    
}

