## Add pi method for lognormal objects
## add back prediction using add_ci_lognormal
add_pi_lm_log <- function(tb, fit, alpha = 0.05, names = NULL) {
 
    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    }
 
  out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
  ##  if(is.null(tb[["pred"]]))
  ##      tb[["pred"]] <- exp(out[, 1])
  tb[[names[1]]] <- exp(out[, 2])
  tb[[names[2]]] <- exp(out[, 3])
  tibble::as_data_frame(tb)
  
}


