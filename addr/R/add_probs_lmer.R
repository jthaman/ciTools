#' Response Probabilities for Linear Mixed Models
#'
#' This function is one of the methods for \code{add_probs}, and is
#' called automatically when \code{add_probs} is used on a \code{fit}
#' of class \code{lmerMod}. It is recommended that one calculate
#' parametric probabilties when modeling with a random intercept
#' mixed model. Otherwise, probabilities may be simulated.
#'
#' @param tb A tibble or Data Frame.
#' @param fit An object of class \code{lmerMod}.
#' @param name NULL or character vector of length two. If \code{NULL},
#'     probabilities will automatically be named by \code{add_pi},
#'     otherwise, the probabilities will be named \code{name} in the
#'     returned tibble.
#' @param q A double. A quantile of the response variable
#' @param type A string, either \code{"parametric"} or \code{"sim"}
#' @param includeRanef A logical. Set whether the predictions and
#'     intervals should be made conditional on the random effects. If
#'     \code{FALSE}, random effects will not be included.
#' @param nSims A positive integer. If \code{type = "sim"}
#'     \code{nSims} will determine the number of simulated draws to
#'     make.
#' @param comparison A character vector of length one. If
#'     \code{comparison = "<"}, then Pr(Y|x < q) is calculated for
#'     each observation in \code{tb}. Must be "<" or ">" for linear,
#'     log-linear and linear mixed models. If \code{fit} is a glm,
#'     then \code{comparison} may also be "<=", ">=", or "=".
#' @param log_response A logical. Set to \code{TRUE} if your model is
#'     a log-linear mixed model.
#' @return A tibble, \code{tb}, with probabilities attached.
#'
#' 
#' @export


add_probs.lmerMod <- function(tb, fit, 
                              q, type = "parametric", 
                              includeRanef = TRUE, name = NULL,
                              comparison = "<", nSims = 200, log_response = FALSE, ...) {
  
    if (is.null(name) && comparison == "<")
        name <- paste("prob_less_than", q, sep="")
    if (is.null(name) && comparison == ">")
        name <- paste("prob_greater_than", q, sep="")

    if (log_response)
        q <- log(q)

    if ((name %in% colnames(tb))) {
        warning ("These Probabilities may have already been appended to your dataframe")
        return(tb)
    }

    if(type == "bootstrap") 
        stop ("this Type is not yet implemented")
    else if (type == "parametric") 
        parametric_probs_mermod(tb, fit, q, name, includeRanef, comparison, ...)
    else if (type == "sim") 
        sim_probs_mermod(tb, fit, q, name, includeRanef, comparison, nSims, ...)
    else  
        stop("Incorrect type specified!")
    
}

parametric_probs_mermod <- function(tb, fit, q, name, includeRanef, comparison){
    
    rdf <- get_resid_df_mermod(fit)
    seGlobal <- get_pi_mermod_var(tb, fit, includeRanef)
    
    if(includeRanef)
        re.form <- NULL
    else
        re.form <- NA

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    
    t_quantile <- (q - tb[["pred"]]) / seGlobal

    if (comparison == "<")
        t_prob <- pt(q = t_quantile, df = rdf)
    if (comparison == ">")
        t_prob <- 1 - pt(q = t_quantile, df = rdf)

    tb[[name]] <- t_prob
    as_data_frame(tb)
}


sim_probs_mermod <- function(tb, fit, q, name, includeRanef, comparison, nSims = 200) {

    if (includeRanef) {
        which <-  "full"
        re.form <- NULL
    } else {
        which <- "fixed"
        re.form <- NA
    }

    pi_out <- predictInterval(fit, tb, which = which, level = 0.95,
                              n.sims = nSims,
                              stat = "median",
                              include.resid.var = TRUE,
                              returnSims = TRUE)
    
    store_sim <- attributes(pi_out)$sim.results
    probs <- apply(store_sim, 1, FUN = calc_prob, quant = q, comparison = comparison)

    if(is.null(tb[["pred"]]))
        tb[["pred"]] <- predict(fit, tb, re.form = re.form)
    tb[[name]] <- probs
    as_data_frame(tb)
    
}
