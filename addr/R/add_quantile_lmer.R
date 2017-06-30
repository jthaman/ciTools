## add quantile method for lmer objects
add_quantile.lmerMod <- function(tb, fit, 
                          prob, includeRanef = TRUE,
                          quantileName = NULL, ...) {

    if (prob <= 0 || prob >= 1)
        stop ("prob should be in (0,1)")
    if (is.null(quantileName)){
        quantileName <- paste("quantile", prob, sep = "")
    }
    if ((quantileName %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    parametric_quantile_mermod(tb, fit, prob, quantileName, includeRanef, ...)
 }

parametric_quantile_mermod <- function(tb, fit, prob, quantileName, includeRanef){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    tb[[quantileName]] <- tb[["pred"]] + qt(prob ,df = rdf) * seGlobal
    tb
}


