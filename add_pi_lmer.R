## add pi method for lmer objects
add_pi.lmerMod <- function(tb, fit, 
                          alpha = 0.05, piType = "parametric", condition_RE = TRUE,
                          includeRanef = TRUE, piNames = NULL, ...) {
    if (is.null(piNames)){
        piNames[1] <- paste("LPB-", alpha/2, sep = "")
        piNames[2] <- paste("UPB-", 1 - alpha/2, sep = "")
    }
    if ((piNames[1] %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe")
        return(tb)
    }

    if(piType == "sim") 
        sim_pi_mermod(tb, fit, alpha, piNames, condition_RE, ...)
    else if(piType == "parametric")
        parametric_pi_mermod(tb, fit, alpha, piNames, includeRanef, condition_RE, ...)
    else
        stop("Incorrect type specified!")

 }

## could make this the default procedure with warnings
parametric_pi_mermod <- function(tb, fit, alpha, piNames, includeRanef, condition_RE){
    if (condition_RE == TRUE)
        reform = NULL
    else
        reform = NA

    X <- model.matrix(reformulate(attributes(terms(fit))$term.labels), tb)
    vcovBetaHat <- vcov(fit)
    
    seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
        diag() %>%
        sqrt()
    
    seRandom <- arm::se.ranef(fit)[[1]][1,]
    rdf <- nrow(model.matrix(fit)) - length(fixef(fit)) -
        (length(attributes(summary(fit)$varcor)$names) + 1)
    se_residual <- sigma(fit)
    if(includeRanef)
        seGlobal <- sqrt(seFixed^2 + seRandom^2 + se_residual^2)
    else
        seGlobal <- sqrt(seFixed^2 + se_residual^2)
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = reform) 
    tb[[piNames[1]]] <- tb[["pred"]] + qt(alpha/2, df = rdf) * seGlobal
    tb[[piNames[2]]] <- tb[["pred"]] + qt(1 - alpha/2, df = rdf) * seGlobal
    tb
    
}

sim_pi_mermod <- function(tb, fit, alpha, piNames, condition_RE, nSims = 1000) {

    if (condition_RE) {
        which = "full"
        reform = NULL
    } else {
        which = "fixed"
        reform = NA
    }
    pi_out <- predictInterval(fit, tb, which = which, level = 1 - alpha,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = TRUE)
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = reform)
    tb[[piNames[1]]] <- pi_out$lwr
    tb[[piNames[2]]] <- pi_out$upr
    tb
    
}
