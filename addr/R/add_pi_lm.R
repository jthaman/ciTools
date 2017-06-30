##add_pi method for lm
add_pi.lm <- function(tb, fit, alpha = 0.05, piNames = NULL, log_response = FALSE){
    
    if (log_response)
        add_pi_lm_log(tb, fit, alpha, piNames)

    else {
        if (is.null(piNames)){
            piNames[1] <- paste("LPB", alpha/2, sep = "")
            piNames[2] <- paste("UPB", 1 - alpha/2, sep = "")
        }
        if ((piNames[1] %in% colnames(tb))) {
            warning ("These PIs may have already been appended to your dataframe")
            return(tb)
        }
        out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
        if(is.null(tb[["pred"]]))
            tb[["pred"]] <- out[, 1]
        if (is.null(tb[[piNames[1]]]))
            tb[[piNames[1]]] <- out[, 2]
        if (is.null(tb[[piNames[2]]]))
            tb[[piNames[2]]] <- out[, 3]
        return(tb)
    }
}

