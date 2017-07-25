get_x_matrix_mermod <- function(tb, fit){
  
  ob <- fit@frame
  cols2add <- names(ob)[!(names(ob) %in% names(tb))]
  for(cn in cols2add){
    tb[[cn]] <- 0
  }
  model.matrix(reformulate(attributes(terms(fit))$term.labels), bind_rows(tb, ob))[1:nrow(tb), ]
  
}


get_prediction_se_mermod <- function(tb, fit){
  
  X <- get_x_matrix_mermod(tb, fit)
  vcovBetaHat <- vcov(fit) %>%
    as.matrix
  X %*% vcovBetaHat %*% t(X) %>% 
    diag %>%
    sqrt
}

make_formula <- function(fixedEffects, randomEffects, rvName = "y"){
    
    fixedPart <- paste(fixedEffects, collapse = "+")
    randomPart <- paste("+ (1|", randomEffects, ")")
    formula(paste(c(rvName, " ~ ", fixedPart, randomPart), collapse = ""))
    
}

add_predictions2 <- function (data, model, var = "pred", ...) {
    data[[var]] <- stats::predict(model, data, ...)
    data
}


get_resid_df_mermod <- function(fit){
    nrow(model.matrix(fit)) - length(fixef(fit)) - 
        (length(attributes(summary(fit)$varcor)$names) + 1)
}

## not useful, consider removal
get_residual_se <- function(fit){
    fit %>%
        VarCorr %>%
        as.data.frame %>%
        last %>%
        last
}


get_pi_mermod_var <- function(tb, fit, includeRanef){
    seFixed <- get_prediction_se_mermod(tb, fit)
    seG <- arm::se.ranef(fit)[[1]][1,]
    sigmaG <- as.data.frame(VarCorr(fit))$sdcor[1]
    se_residual <- sigma(fit)
    
    if(includeRanef)
        return(sqrt(seFixed^2 + seG^2 + se_residual^2))
    else
        return(sqrt(seFixed^2 + sigmaG^2 + se_residual^2))
}


calc_prob <- function(x, quant, comparison){
    if (comparison == "<")
        mean(x < quant)
    else if (comparison == ">")
        mean(x > quant)
    else if (comparison == "<=")
        mean(x <= quant)
    else if (comparison == ">=")
        mean(x >= quant)
    else if (comparison == "=")
        mean (x == quant)
    else
        stop ("Malformed probability statement, comparison must be <, >, =, <=, or >=")
}

my_pred_full <- function(fit) {
    predict(fit, newdata=.tb_temp1234567890, re.form = NULL)
}

my_pred_fixed <- function(fit) {
    predict(fit, newdata=.tb_temp1234567890, re.form = NA)
}

boot_quants <- function(merBoot, alpha) {
    return(
        data.frame(fit = apply(merBoot$t, 2, quantile, probs = 0.5),
                   lwr = apply(merBoot$t, 2, quantile, probs = alpha / 2),
                   upr = apply(merBoot$t, 2, quantile, probs = 1 - alpha / 2)))
}

