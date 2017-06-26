## add_ci method for GLMs
## TODO : add Bootstrap method (not a primary concern)
## TODO : Add more families

add_ci.glm <- function(tb, fit, alpha = 0.05, ciNames = NULL,
                       type = "response", ciType = "parametric", bootReplicates = 1000){

    if (is.null(ciNames)){
        ciNames[1] <- paste("LCB-", alpha/2, sep = "")
        ciNames[2] <- paste("UCB-", 1 - alpha/2, sep = "")
    }
    if ((ciNames[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe")
        return(tb)
    }

    if ( !(fit$fit$family %in% c("normal", "poisson", "gamma", "binomial")))
        stop ("We don't support this glm family'")

    if (ciType == "bootstrap")
        stop ("not yet implemented")
    else if (ciType == "parametric")
        parametric_ci_glm(tb, fit, alpha, ciNames, type,)
    else
        if(!(ciType %in% c("bootstrap", "parametric"))) stop("Incorrect interval type specified!")

}

parametric_ci_glm <- function(tb, fit, alpha, ciNames, type = "response"){

    inverselink <- fit$family$linkinv
    if(type == "response"){
        out <- glm_ci_response(tb, fit, alpha, inverselink) %>%
            plyr::rename(c("lwr" = ciNames[1], "upr" = ciNames[2])) %>%
            select(-.se)
    }

    if(type == "link"){
        out <- glm_ci_linear_predictor(tb, fit, alpha) %>%
            plyr::rename(c("lwr" = ciNames[1], "upr" = ciNames[2])) %>%
            select(-.se)
    }
    if(!(type %in% c("response", "link"))) top("Incorrect interval type specified!")
    out
}

glm_ci_response <- function(tb, fit, alpha, ilink){

    if(is.null(tb[["pred"]]))
        tb <- modelr::add_predictions(tb, fit)
    tb %>%
        add_standard_error(fit) %>%
        mutate(
            lwr = ilink(pred + qt(alpha/2, df = fit$df.residual) * .se),
            upr = ilink(pred + qt(1 - alpha/2, df = fit$df.residual) * .se)
        )

}

glm_ci_linear_predictor <- function(tb, fit, alpha){

    if(is.null(tb[["pred"]]))
        tb <- modelr::add_predictions(tb, fit)
    tb %>%
        add_standard_error(fit) %>%
        mutate(
            lwr = pred + qt(alpha/2, df = fit$df.residual) * .se,
            upr = pred + qt(1 - alpha/2, df = fit$df.residual) * .se
        )

}
