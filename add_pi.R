# Add PI plus generics

# Defining the generic for "add_pi"
add_pi <- function(tb, fit, alpha = 0.05, piNames = c("LPB", "UPB"), ...){
  UseMethod("add_pi")
}

## add pi method for GLM

## add_pi method for LogNormal Models

## add_pi method for GLMM
add_pi.lmerMod <- function(tb, fit2, 
                           alpha = 0.05, piType = "parametric", 
                           includeRanef = T, piNames = c("LCB", "UCB"), ...){
  if(piType == "bootstrap") return(
    bootstrap_pi_mermod(tb, fit2, alpha, piNames, ...)
  ) else
    if(piType == "parametric") return(
      parametric_pi_mermod(tb, fit2, alpha, piNames, includeRanef)
    ) else
      if(!(piType %in% c("bootstrap", "parametric"))) stop("Incorrect type specified!")
  
}

parametric_pi_mermod <- function(tb, fit, alpha, piNames, includeRanef){
    X <- model.matrix(fit)
    vcovBetaHat <- vcov(fit)
    
    seFixed <- X %*% vcovBetaHat %*% t(X) %>% 
        diag() %>%
        sqrt()
    
    seRandom <- arm::se.ranef(fit)[[1]][1,]
    rdf <- nrow(model.matrix(fit)) - length(fixef(fit)) - (length(attributes(summary(fit)$varcor)$names) + 1)

    seResidual <- fit %>%
        VarCorr %>%
        as.data.frame %>%
        last %>%
        last

    if(includeRanef){seGlobal <- sqrt(seFixed^2 + seRandom^2 + seResidual^2)
        }else{seGlobal <- seFixed}
    
    tb[[piNames[1]]] <- tb[["pred"]] + qt(alpha/2,df = rdf) * seGlobal
    if(is.null(tb[["pred"]])) tb <- modelr::add_predictions(tb, fit)
    tb[[piNames[2]]] <- tb[["pred"]] + qt(1 - alpha/2, df = rdf) * seGlobal
    tb
}

#######testing#########
fm1 <- lmer(
    formula = distance ~ age*Sex + (age|Subject)
    , data = Orthodont
) 

comp.data <- parametric_pi_mermod(tb= Orthodont, fit = fm1, alpha = 0.05, piNames = c("LPB", "UPB"), includeRanef = TRUE)
comp.data <- select(comp.data, c(pred, LPB, UPB))
comp.data <- comp.data[1:20,]

comp.data.new <- predictInterval(fm1, newdata = Orthodont, include.resid.var = FALSE)
comp.data.new <- comp.data.new[c(1,3,2)]
comp.data.new <- comp.data.new[1:20,]

combined <- cbind(comp.data, comp.data.new)

p <- ggplot(combined) + geom_line(aes(x=1:20, y=pred, colour="red")) 
p + geom_ribbon(data=combined,aes(x=1:20, y=pred, ymin=LPB,ymax=UPB),alpha=0.3) + 
geom_ribbon(data=combined,aes(x=1:20, y=pred, ymin=lwr,ymax=upr),alpha=0.1)


data(Orthodont, package="nlme")
fm1 <- lmer(distance ~ age + (age|Subject), data = Orthodont)
(vc <- VarCorr(fm1))  ## default print method: standard dev and corr
## both variance and std.dev.
print(vc,comp=c("Variance","Std.Dev."),digits=2)
## variance only
print(vc,comp=c("Variance"))
as.data.frame(vc)
as.data.frame(vc,order="lower.tri")
