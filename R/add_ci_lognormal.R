# Copyright (C) 2017 Institute for Defense Analyses
#
# This file is part of ciTools.
#
# ciTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ciTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ciTools. If not, see <http://www.gnu.org/licenses/>.

## Special thanks to Dr. Rebecca Dickins for her derivations of the
## confidence intervals of E[Y|x] in a log-Normal model.

get_sigma_mle <- function(tb, fit, alpha){
    X <- model.matrix(fit)
    n <- NROW(X)
    out <- predict(fit, tb, se.fit = TRUE, interval = "confidence", level = 1 - alpha)
    rmse <- out$residual.scale
    residual_df <- out$df
    sigma_mle <- sqrt(rmse^2 * residual_df / n)
    sigma_mle
}

create_cov_mat <- function(sigma_mle, tb, fit){
    X <- model.matrix(fit)
    n <- NROW(X)
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
    return(list(pred = pred, Xpred = Xpred, Xbeta = Xbeta))
}


get_se_pred <- function(pred, Xpred, V, p){
    se_pred <- rep(0, length(pred))
    ## this would be worth rewriting
    for(i in 1:length(pred)){
        derivative_vector <- c(cbind(Xpred[i,] * pred[i]), cbind(qnorm(p) * pred[i]))
        se_pred[i] <- sqrt(t(derivative_vector) %*% V %*% derivative_vector)
    }
    se_pred
}

add_ci_lm_log <- function(tb, fit, alpha = 0.05, names = NULL, yhatName){

    if (is.null(names)) {
        names[1] <- paste("LCB", alpha/2, sep = "")
        names[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe. Overwriting.")
    }
    sigma_mle <- get_sigma_mle(tb, fit, alpha)
    V <- create_cov_mat(sigma_mle, tb, fit)
    p <- pnorm(sigma_mle / 2)
    pred_list <- get_Xpred_list(tb, fit, p, sigma_mle)
    pred <- pred_list$pred
    Xpred <- pred_list$Xpred
    Xbeta <- pred_list$Xbeta
    z_quantile <- qnorm(1 - alpha/2) 
    se_pred <- get_se_pred(pred, Xpred, V, p)

    lwr <- exp(Xbeta + qnorm(p) * sigma_mle - z_quantile * se_pred / pred)
    upr <- exp(Xbeta + qnorm(p) * sigma_mle + z_quantile * se_pred / pred)
    
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- pred
    tb[[names[1]]] <- lwr
    tb[[names[2]]] <- upr

    tibble::as_data_frame(tb)
}


