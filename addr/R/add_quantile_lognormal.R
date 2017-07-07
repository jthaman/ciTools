## add quantile function for log-linear models
## add back predictions using add_ci_lognormal
add_quantile_lm_log <- function(tb, fit, p, name = NULL) {
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)) {
        name <- paste("quantile", p, sep = "")
    }

    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    
    out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
    fitted <- out$fit[,1]
    residual_df <- out$df
    se_fitted <- out$se.fit
    resid_var <- out$residual.scale^2
    se_pred <- sqrt(resid_var + se_fitted^2)
    t_quantile <- qt(p = p, df = residual_df)
    out_quantiles <- exp(fitted + se_pred * t_quantile)

    #if(is.null(tb[["pred"]]))
       # tb[["pred"]] <- exp(fitted)
    tb[[name]] <- out_quantiles
    as_data_frame(tb)
}


