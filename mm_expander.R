library(tidyverse)

## This only really works well with very clean data
mm_expander <- function(tb, fit){
    lst <- list()
    for (i in colnames(tb)){
        if (is.double(tb[[i]]))
            lst[[i]] <- seq(min(tb[[i]]), max(tb[[i]]), length.out = 101)
        else if (is.character(tb[[i]]))
            lst[[i]] <- unique(tb[[i]])
        else if (is.integer(tb[[i]]))
            lst[[i]] <- seq(min(tb[[i]]), max(tb[[i]]), by = 1)
        else if (is.logical(tb[[i]]))
            lst[[i]] <- c(TRUE, FALSE)
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

