#' @export

add_probs.lmerMod <- function(tb, fit, 
                              quant, probType = "parametric", 
                              includeRanef = TRUE, probName = NULL,
                              comparison = "<", nSims = 200, log_response = FALSE, ...) {
    if (log_response)
        quant <- log(quant)

    if (is.null(probName) && comparison == "<")
        probName <- paste("PrYless", quant, sep="")
    if (is.null(probName) && comparison == ">")
        probName <- paste("PrYgreater", quant, sep="")

    if ((probName %in% colnames(tb))) {
        warning ("These Probabilities may have already been appended to your dataframe")
        return(tb)
    }

    if(probType == "bootstrap") 
        stop ("this Type is not yet implemented")
    else if (probType == "parametric") 
        parametric_probs_mermod(tb, fit, quant, probName, includeRanef, comparison, ...)
    else if (probType == "sim") 
        sim_probs_mermod(tb, fit, quant, probName, includeRanef, comparison, nSims, ...)
    else  
        stop("Incorrect type specified!")
    
}

parametric_probs_mermod <- function(tb, fit, quant, probName, includeRanef, comparison){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    
    t_quantile <- (quant - tb[["pred"]]) / seGlobal

    if (comparison == "<")
        t_prob <- pt(q = t_quantile, df = rdf)
    if (comparison == ">")
        t_prob <- 1 - pt(q = t_quantile, df = rdf)

    tb[[probName]] <- t_prob
    as_data_frame(tb)
}


sim_probs_mermod <- function(tb, fit, quant, probName, includeRanef, comparison, nSims = 200) {

    if (includeRanef) {
        which <-  "full"
        re.form <- NULL
    } else {
        which <- "fixed"
        re.form <- NA
    }

    pi_out <- predictInterval(fit, tb, which = which, level = 0.95,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = TRUE,
                              returnSims = TRUE)
    
    store_sim <- attributes(pi_out)$sim.results
    probs <- apply(store_sim, 1, FUN = calc_prob, quant = quant, comparison = comparison)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    tb[[probName]] <- probs
    as_data_frame(tb)
    
}
