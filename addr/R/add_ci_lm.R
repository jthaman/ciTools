#add_ci method for lm
add_ci.lm <- function(tb, fit, alpha = 0.05, ciNames = NULL){
    
    if (is.null(ciNames)){
        ciNames[1] <- paste("LCB", alpha/2, sep = "")
        ciNames[2] <- paste("UCB", 1 - alpha/2, sep = "")
    }
    if ((ciNames[1] %in% colnames(tb))) {
        warning ("These CIs may have already been appended to your dataframe")
        return(tb)
    }
    out <- predict(fit, tb, interval = "confidence", level = 1 - alpha)
    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- out[, 1]
    if (is.null(tb[[ciNames[1]]]))
        tb[[ciNames[1]]] <- out[, 2]
    if (is.null(tb[[ciNames[2]]]))
        tb[[ciNames[2]]] <- out[, 3]
    return(tb)
    
}

