library(lme4)
library(arm)
library(merTools)

J <- 10
n <- J*(J+1)/2
group <- rep (1:J, 1:J)
mu.a <- 5
sigma.a <- 2
a <- rnorm (J, mu.a, sigma.a)
b <- -3
x <- rnorm (n, 2, 1)
sigma.y <- 6
y <- rnorm (n, a[group] + b*x, sigma.y)
u <- runif (J, 0, 3)
tb <- as_data_frame(cbind (y, x, group))

mm <- lmer (y ~ x + (1 | group))
u.full <- u[group]
mm <- lmer (y ~ x + u.full + (1 | group)

##tb <- sleepstudy

##fit <- mm <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)
##mm <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)

Ran_True <- add_ci.lmerMod(tb, mm, includeRanef = TRUE, condition_RE = TRUE)

Ran_False <- add_ci.lmerMod(tb, mm, includeRanef = FALSE, condition_RE = FALSE)

Ran_True_sim <- add_ci.lmerMod(tb, mm, alpha = 0.05, ciType = "sim", condition_RE = TRUE, nSims = 1000)


Ran_False_sim <- add_ci.lmerMod(tb, mm, alpha = 0.05, ciType = "sim", condition_RE = FALSE, nSims = 1000)

comp.data <- rbind(
    data.frame(Predict.Method="wald_w_ran", x=(1:nrow(Ran_True))-0.1, Ran_True),
    data.frame(Predict.Method="wald_wo_ran", x=(1:nrow(Ran_False))-0.2, Ran_False),
    data.frame(Predict.Method="PI_w_ran", x=(1:nrow(Ran_True_sim))+0.1, Ran_True_sim),
    data.frame(Predict.Method="PI_wo_ran", x=(1:nrow(Ran_False_sim))+0.2, Ran_False_sim)
)

ggplot(aes(x=x, y=pred, ymin=LCB.0.025, ymax=UCB.0.975, color=Predict.Method), data=comp.data[c(1:30,
                                                                                    181:210,
                                                                                    181:210+180,
                                                                                    181:210+360
                                                                                    ),]) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
  scale_color_brewer(type = "qual", palette = 2)



ggplot(aes(x=x, y=pred, ymin=LCB.0.025, ymax=UCB.0.975, color=Predict.Method), data=comp.data) +
  geom_point() + 
  geom_linerange() +
  labs(x="Index", y="Prediction w/ 95% PI") +
  theme_bw() +  theme(legend.position="bottom") +
scale_color_brewer(type = "qual", palette = 2) +
geom_point(aes(x = x, y = y), color = "black")
