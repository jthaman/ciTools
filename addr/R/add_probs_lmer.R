#' @export

add_probs.lmerMod <- function(tb, fit, 
                              quant, type = "parametric", 
                              includeRanef = TRUE, name = NULL,
                              comparison = "<", nSims = 200, log_response = FALSE, ...) {
    if (log_response)
        quant <- log(quant)

    if (is.null(name) && comparison == "<")
        name <- paste("Pr_less_", quant, sep="")
    if (is.null(name) && comparison == ">")
        name <- paste("Pr_greater_", quant, sep="")

    if ((name %in% colnames(tb))) {
        warning ("These Probabilities may have already been appended to your dataframe")
        return(tb)
    }

    if(type == "bootstrap") 
        stop ("this Type is not yet implemented")
    else if (type == "parametric") 
        parametric_probs_mermod(tb, fit, quant, name, includeRanef, comparison, ...)
    else if (type == "sim") 
        sim_probs_mermod(tb, fit, quant, name, includeRanef, comparison, nSims, ...)
    else  
        stop("Incorrect type specified!")
    
}

parametric_probs_mermod <- function(tb, fit, quant, name, includeRanef, comparison){
    
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

    tb[[name]] <- t_prob
    as_data_frame(tb)
}


sim_probs_mermod <- function(tb, fit, quant, name, includeRanef, comparison, nSims = 200) {

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
    tb[[name]] <- probs
    as_data_frame(tb)
    
}
