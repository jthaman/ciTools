tb <- sleepstudy
mm <- lmer(Reaction ~ Days + (1 | Subject), data = sleepstudy)

simulation_ci_mermod(tb, mm, 0.1, ciNames = c("LCB", "UCB"), includeRanef=FALSE, nSims = 1000)

simulation_ci_mermod(tb, mm, 0.1, ciNames = c("LCB", "UCB"), includeRanef=TRUE, nSims = 1000)
