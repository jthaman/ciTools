#add_ci method for GLMs
## TODO : add Matt's edits
## TODO : add Bootstrap method (not a primary concern)

add_ci.glm <- function(tb, fit, alpha = 0.05, ciNames = c("LCB", "UCB"),
                       type = "response", ciType = "parametric", bootReplicates = 1000){

  if(ciType == "bootstrap") return(
    bootstrap_ci_glm(tb, fit, alpha, ciNames, type, bootReplicates, ...)
  ) else
    if(ciType == "parametric") return(
      parametric_ci_glm(tb, fit, alpha, ciNames, type,)
    ) else
      if(!(ciType %in% c("bootstrap", "parametric"))) stop("Incorrect interval type specified!")

}

## add_ci Wald method for glm

parametric_ci_glm <- function(tb, fit, alpha, ciNames, type, ...){
    out <- predict(fit, tb, se.fit = TRUE)
    crit_val <- qt(p = 1 - alpha/2, df = fit$df.residual)
    inverselink <- fit$family$linkinv
    if(type == "response"){
        upr <- inverselink(out$fit + crit_val * out$se.fit)
        lwr <- inverselink(out$fit - crit_val * out$se.fit)
        pred <- inverselink(out$fit)
    }else{
        upr <- out$fit + crit_val * out$se.fit
        lwr <- out$fit - crit_val * out$se.fit
        pred <- out$fit
    }
    if(is.null(tb[["pred"]])) tb[["pred"]] <- pred
    tb[[ciNames[1]]] <- lwr
    tb[[ciNames[2]]] <- upr
    tb

}

