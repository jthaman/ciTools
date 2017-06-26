# Draw figure 1 from Escobar and Meeker, 1992.
library(survival)

fit <- survreg(Surv(time,status) ~ age + I(age^2), data=stanford2, 
	dist='lognormal') 
pred <- predict(fit, newdata=list(age=1:65), se.fit = TRUE, type='response') 

# Predicted Weibull survival curve for a lung cancer subject with
#  ECOG score of 2
lfit <- survreg(Surv(time, status) ~ ph.ecog, data=lung)
pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(lfit, newdata=data.frame(ph.ecog=2), type='quantile',
                 p=pct, se=TRUE)
matplot(cbind(ptime$fit, ptime$fit + 2*ptime$se.fit,
                         ptime$fit - 2*ptime$se.fit)/30.5, 1-pct,
        xlab="Months", ylab="Survival", type='l', lty=c(1,2,2), col=1)
