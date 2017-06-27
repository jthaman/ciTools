## add_ci method for lme4/merMod objects

add_ci.lmerMod <- function(tb, fit, 
                           alpha = 0.05, ciType = "parametric", condition_RE = TRUE,
                           includeRanef = TRUE, ciNames = NULL, nSims = 1000, ...){

    if (is.null(ciNames)){
        ciNames[1] <- paste("LCB-", alpha/2, sep = "")
        ciNames[2] <- paste("UCB-", 1 - alpha/2, sep = "")
    }
    if ((ciNames[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe")
        return(tb)
    }

    if (ciType == "bootstrap") 
        bootstrap_ci_mermod(tb, fit, alpha, ciNames, condition_RE, nSims, ...)
    else if (ciType == "parametric")
        parametric_ci_mermod(tb, fit, alpha, ciNames, includeRanef, condition_RE, ...)
    else if (ciType == "sim")
        sim_ci_mermod(tb, fit, alpha, ciNames, condition_RE, nSims, ...)

    else
        stop("Incorrect type specified!")
    
}

## this function should be used with exterme caution
parametric_ci_mermod <- function(tb, fit, alpha, ciNames, includeRanef, condition_RE){
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
    if(includeRanef)
        seGlobal <- sqrt(seFixed^2 + seRandom^2)
    else
        seGlobal <- seFixed
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = reform) 
    tb[[ciNames[1]]] <- tb[["pred"]] + qt(alpha/2, df = rdf) * seGlobal
    tb[[ciNames[2]]] <- tb[["pred"]] + qt(1 - alpha/2, df = rdf) * seGlobal
    tb
}

## simulation method using predictInterval
sim_ci_mermod <- function(tb, fit, alpha, ciNames, condition_RE, nSims = 1000) {

    if (condition_RE) {
        which = "full"
        reform = NULL
    } else {
        which = "fixed"
        reform = NA
    }

    ci_out <- predictInterval(fit, tb, which = which, level = 1 - alpha,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = FALSE)

    if(is.null(tb[["pred"]])) 
        tb[["pred"]] <- predict(fit, tb, re.form = reform) 
    tb[[ciNames[1]]] <- ci_out$lwr
    tb[[ciNames[2]]] <- ci_out$upr
    tb
    
}

mySumm <- function(.) {
    predict(., newdata=tb, re.form = NULL)
}

sumBoot <- function(merBoot, alpha = alpha) {
    return(
        data.frame(fit = apply(merBoot$t, 2, quantile, probs = 0.5),
                   lwr = apply(merBoot$t, 2, quantile, probs = alpha / 2),
                   upr = apply(merBoot$t, 2, quantile, probs = 1 - alpha / 2)
                   )
    )
}


## only conditional estimation is supported right now
bootstrap_ci_mermod <- function(tb, fit, alpha, ciNames, condition_RE = TRUE, nSims) {
    
    if (condition_RE) { 
        rform = NULL
    } else {
        rform = NA
    }
        
    boot_obj <- lme4::bootMer(fit, mySumm, nsim=nSims, use.u=FALSE, type="parametric", rform)

    ci_out <- sumBoot(boot_obj, alpha) 

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- ci_out$fit
    tb[[ciNames[1]]] <- ci_out$lwr
    tb[[ciNames[2]]] <- ci_out$upr
    tb
    
}

## ##Bootstrap CIs for merMod objects
## ##Note that the bootstrapping is done on the random effects only. 

## oldbootstrap_ci_mermod <- function(tb, fit, alpha, ciNames, nBS = 999){
    
##     X <- model.matrix(fit)
##     vcovBetaHat <- vcov(fit)
    
##     seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
##         diag() %>%
##         sqrt()
    
##     wholePlotBS <- sample_ranefs(fit, nBS, alpha) 
##     if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
##     tb[[ciNames[1]]] <- tb[["pred"]] + qnorm(alpha/2) * seFixed + wholePlotBS[1]
##     tb[[ciNames[2]]] <- tb[["pred"]] + qnorm(1 - alpha/2) * seFixed + wholePlotBS[2]
##     tb
    
## }

## sample_ranefs <- function(fit, nBS, alpha){
##     if(length(ranef(fit)[[1]][[1]]) < 5) warning(
##                                              "Fewer than 5 whole plots; consider using parametric methods")
##     ranef(fit)[[1]][[1]] %>% 
##         base::sample(size = nBS, replace = T) %>%
##         quantile(probs = c(alpha/2, 1-alpha/2))
## }
