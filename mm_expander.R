## This only really works well with very clean data
## This function presently does not incorporate anything about the fit

mm_expander <- function(tb, fit){
    tb <- na.omit(tb) ## is this okay?
    pred_orig <- predict(fit, expanded)
    tb <- mutate(tb, pred_orig)
    lst <- list()
    for (i in colnames(tb)){
        if (is.factor(tb[[i]]))
            lst[[i]] <- unique(tb[[i]])
        else if (is.double(tb[[i]]))
            lst[[i]] <- seq(min(tb[[i]]), max(tb[[i]]), length.out = 101)
        else if (is.character(tb[[i]]))
            lst[[i]] <- unique(tb[[i]])
        else if (is.integer(tb[[i]]))
            lst[[i]] <- seq(min(tb[[i]]), max(tb[[i]]), by = 1)
        else {
            next
            warning("There is a missing value or a data type that I cannot handle")
        }
    }
    expanded <- lst %>% crossing_
    pred <- predict(fit, expanded)
    expanded_pred <- mutate(expanded, pred)
    rbind(tb, expanded_pred)
}

