## add ci for log linear models
## add option to use this since it's not a method

get_sigma_mle <- function(tb, fit){
    X <- model.matrix(fit)
    out <- predict(fit, tb, se.fit = TRUE, interval = "confidence", level = 1 - alpha)
    rmse <- out$residual.scale
    residual_df <- out$df
    sigma_mle <- sqrt(rmse^2 * error.df / n)
    sigma_mle
}

create_cov_mat <- function(sigma_mle, tb, fit){
    X <- model.matrix(fit)
    sigma_var_mle <- sigma_mle^2 /(2 * n)
    beta <- fit$coefficients
    beta_cov_mat <- sigma_mle^2 * solve(t(X) %*% X) 
    k <- length(beta)
    V <- matrix(0, nrow = k + 1, ncol = k + 1 )
    V[1:k,1:k] <- beta_cov_mat
    V[k+1,k+1] <- sigma_var_mle
    V
}

get_Xpred_list <- function(tb, fit, p, sigma_mle){
    Xbeta <- predict(fit, tb)
    pred <- exp(Xbeta + qnorm(p) * sigma_mle)
    reg_mat <- cbind(tb, pred)
    names(reg_mat)[ncol(reg_mat)] = "Prediction"
    lhs <- "Prediction"
    rhs <- as.character(fit$call$formula)[3]
    form <- paste(lhs, "~", rhs)
    fitPred <- lm(formula = form, data = reg_mat)
    Xpred <- model.matrix(fitPred)
    return(list(pred = pred, Xpred = Xpred))
}


get_se_pred <- function(pred, Xpred, V, p){
    se_pred <- rep(0, length(pred))
    for(i in 1:length(pred)){
        derivative_vector <- c(cbind(Xpred[i,] * pred[i]), cbind(qnorm(p) * pred[i]))
        se_pred[i] <- sqrt(t(derivative_vector) %*% V %*% derivative_vector)
    }
    se_pred
}

add_ci_lm_log <- function(tb, fit, alpha = 0.05, ciNames = c("LCB", "UCB")){
    sigma_mle <- get_sigma_mle(tb, fit)
    V <- create_cov_mat(sigma_mle, tb, fit)
    p <- pnorm(sigma_mle / 2)
    pred_list <- get_Xpred_list(tb, fit, p, sigma_mle)
    pred <- pred_list$pred
    Xpred <- pred_list$Xpred
    z_quantile <- qnorm(1 - alpha/2) 
    se_pred <- get_se_pred(pred, Xpred, V, p)

    lwr <- exp(Xbeta + qnorm(p) * sigma_mle - z_quantile * se_pred / pred)
    upr <- exp(Xbeta + qnorm(p) * sigma_mle + z_quantile * se_pred / pred)
    
    if(is.null(tb[["pred"]])) tb[["pred"]] <- pred
    tb[[ciNames[1]]] <- lwr
    tb[[ciNames[2]]] <- upr
    tb
    
}


