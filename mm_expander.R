## This only really works well with very clean data
## TODO : only expand over predictors!
mm_expander <- function(tb, fit, dbls){
    tb <- na.omit(tb) 
    pred <- predict(fit, tb)
    tb <- mutate(tb, pred)
    lst <- list()
    cands <- intersect(colnames(tb), variable.names(fit))
    for (i in cands){
        if (is.factor(tb[[i]]))
            lst[[i]] <- unique(tb[[i]])
        else if (is.double(tb[[i]]))
            lst[[i]] <- seq(min(tb[[i]]), max(tb[[i]]), length.out = dbls)
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
    ##rbind(tb, expanded_pred)
    expanded_pred
}

