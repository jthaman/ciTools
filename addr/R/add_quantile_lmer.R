#' @export

add_quantile.lmerMod <- function(tb, fit, 
                          prob, includeRanef = TRUE,
                          name = NULL, ...) {

    if (prob <= 0 || prob >= 1)
        stop ("prob should be in (0,1)")
    if (is.null(name)){
        name <- paste("quantile", prob, sep = "")
    }
    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe")
        return(tb)
    }
    parametric_quantile_mermod(tb, fit, prob, name, includeRanef, ...)
 }

parametric_quantile_mermod <- function(tb, fit, prob, name, includeRanef){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    tb[[name]] <- tb[["pred"]] + qt(prob ,df = rdf) * seGlobal
    as_data_frame(tb)
}


