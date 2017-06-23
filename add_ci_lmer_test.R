library(lme4)
library(arm)
library(merTools)

tb <- sleepstudy

fit <- mm <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
##mm <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)

Ran_True <- add_ci.lmerMod(tb, mm, includeRanef = TRUE, condition_RE = TRUE)

Ran_False <- add_ci.lmerMod(tb, mm, includeRanef = FALSE, condition_RE = FALSE)

PI_Ran_False <- simulation_ci_mermod(tb, mm, 0.1, ciNames = c("LCB", "UCB"), includeRanef=FALSE, nSims = 1000)

PI_Ran_True <- simulation_ci_mermod(tb, mm, 0.1, ciNames = c("LCB", "UCB"), includeRanef=TRUE, nSims = 1000)

comp.data <- rbind(
    data.frame(Predict.Method="wald_w_ran", x=(1:nrow(Ran_True))-0.1, Ran_True),
    data.frame(Predict.Method="wald_wo_ran", x=(1:nrow(Ran_False))-0.2, Ran_False),
    data.frame(Predict.Method="PI_wo_ran", x=(1:nrow(PI_Ran_False))+0.1, PI_Ran_False),
    data.frame(Predict.Method="PI_w_ran", x=(1:nrow(PI_Ran_True))+0.2, PI_Ran_True)
)

ggplot(aes(x=x, y=pred, ymin=LCB, ymax=UCB, color=Predict.Method), data=comp.data[c(1:30,
                                                                                    181:210,
                                                                                    181:210+180,
                                                                                    181:210+360
                                                                                    ),]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)
