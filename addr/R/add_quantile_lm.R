#' @export

add_quantile.lm <- function(tb, fit, prob, quantileName = NULL, log_response = FALSE){
    if (log_response)
        add_quantile_lm_log(tb, fit, prob, quantileName)
    else {
        if (prob == 0.5)
            warning ("The 0.5 quantile is equal to the fitted values")
        if (prob <= 0 || prob >= 1)
            stop ("prob should be in (0,1)")
        if (is.null(quantileName))
            quantileName <- paste("quantile", prob, sep="")
        if (quantileName %in% colnames(tb)) {
            warning ("These quantiles may have already been appended to your dataframe")
            return(tb)
        }
        out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
        fitted <- out$fit[,1]
        residual_df <- out$df
        se_fitted <- out$se.fit
        resid_var <- out$residual.scale^2
        se_pred <- sqrt(resid_var + se_fitted^2)
        t_quantile <- qt(p = prob, df = residual_df)
        out_quantiles <- fitted + se_pred * t_quantile
        if (is.null(tb[["pred"]]))
            tb[["pred"]] <- fitted
        tb[[quantileName]] <- out_quantiles
        as_data_frame(tb)
    }
}

