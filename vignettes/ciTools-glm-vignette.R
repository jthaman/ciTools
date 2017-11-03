#' ---
#' title: Using ciTools with Generalized Linear Models
#' author: John Haman
#' date: November 2, 2017
#' output:
#'    html_document:
#'      toc: true
#' ---

#+ echo=FALSE
library(tidyverse)
library(ciTools)
library(MASS)
library(arm)
set.seed(20171102)

#' In this vignette we will discuss the current ability of `ciTools`
#' to handle generalized linear models. Small simulations will be
#' provided in addition to examples that show how to use `ciTools` to
#' quantify uncertainty in a GLM. Primarily we focus on the
#' Logistic and Poisson models, but `ciTools`'s method for handling GLMs
#' is not limited to these models. Note that the *Logistic-Binomial*
#' model is handled in a separate vignette. 
#' 
#' ## The Generalized Linear Model
#'
#' Generalized linear models are an extension of linear models that
#' seek to accomodate certain types of non-linear relationships. The
#' manner in which the non-linear relationship is addressed also
#' allows users to perform inferences on data that are not strictly
#' continuous. GLMs are the most common model type that allow for a
#' non-linear relationship between the response variable $y$ and
#' covariates $X$. Recall that linear regression directly predicts a
#' continuous response from the linear predictor $X \beta$. A GLM
#' extends this linear prediction scheme through the following
#' components:
#'
#' 1. A response vector $y = (y_1, \ldots, y_n)$.
#' 
#' 2. A linear predictor $X \beta$
#' 
#' 3. A monotonic and everywhere differentiable link function $g$,
#' which transforms the linear predictor: $\hat{y} = g^{-1}(X \beta)$.
#'
#' 4. A response distribution: $f(y|\mu)$ from the exponential family
#' with mean $\mu = g^{-1} (X \beta)$ and a variance $\mathrm{Var}(y) =
#' f(\mu)$.
#' 
#' Other components, such as over dispersion parameters and off-set
#' terms are possible, but not common to all GLMs. The most common
#' GLMs in practice are the Logistic model (Bernoulli response with
#' logit link) and the Poisson model with log-link. These are detailed
#' below for convenience.
#' 
#' 1. Logistic Regression:
#' $$
#' \begin{equation}
#' \begin{split}
#' & y|X  \sim \mathrm{Binomial}(1, p) \\
#' & g(p) = \log \left( \frac{p}{1-p} \right) = X\beta \\
#' & \mathrm{E}[y|X] = p = \frac{\exp(X \beta)}{ 1 + \exp(X \beta)}
#' \end{split}
#' \end{equation}
#' $$
#' 
#' 2. Poisson Regression with the $\log$ link function:
#' $$
#' \begin{equation}
#' \begin{split}
#' & y|X  \sim \mathrm{Poisson}(\lambda) \\
#' & g(\lambda) = \log \left( \lambda \right) = X\beta \\
#' & \mathrm{E}[y|X] = \lambda = \exp(X \beta)
#' \end{split}
#' \end{equation}
#' $$
#'
#' Due to the variety of options available, fitting generalized linear
#' models is more complicated than fitting linear models. In **R**,
#' `glm` is the starting point for handling GLM fits, and is the
#' currently the only model fitting function that is supported by
#' `ciTools`. We can use `glm` to fit Logistic, Poisson, Quasipoisson,
#' Gamma and Guassian GLMs. All of these models are handled by
#' `ciTools`.
#'
#' # Overview of `ciTools` methods for GLMs
#'
#' Unlike linear models, interval estimates pertaining to GLMs
#' generally do not have clean parametric forms. Parametric interval
#' estimates are available in certain cases, and where ever available,
#' `ciTools` will choose to implement them by default. Below we detail
#' precisely what `ciTools` is doing under the hood when one of the
#' core functions (`add_ci`, `add_pi`, `add_probs`, `add_quantile`) is
#' called on a GLM fit.
#'
#' ## Confidence Intervals
#' 
#' For any model fit by `glm`, `add_ci()` may compute confidence
#' intervals for predictions using either a parametric method or a
#' bootstrap. The parametric confidence interval method computes
#' confidence intervals of the scale of the linear predictor $X \beta$
#' and transforms the intervals to the response level through the
#' inverse link function $g^{-1}$. Confidence intervals on the linear
#' predictor level are computed using a Normal distribution for
#' Logistic and Poisson regressions and a $t$ distribution otherwise.
#'
#' The default method is parametric and is called with `add_ci(data,
#' fit, ...)`. This is the method we generally recommend for
#' constructing confidence intervals for model prediction. The
#' bootstrap method is called with `add_ci(data, fit, type = "boot",
#' ...)` and is included mostly for making comparisons against the
#' parametric method. Note that there are multiple methods of
#' bootstrap for regression models (resampling cases, resampling
#' residuals, parametric, etc.). The bootstrap method employed by
#' `ciTools` in `add_ci.glm()` resamples cases and iteratively refits
#' the model to determine confidence intervals.
#'
#' For comparison, we show an example of the confidence intervals for
#' the probability estimates of a Logistic regression model.
#'
#' 
x <- rnorm(100, mean = 5)
y <- rbinom(n = 100, size = 1, prob = invlogit(-20 + 4*x))
df <- data.frame(x = x, y = y)
fit <- glm(y ~ x, family = binomial)

#' We use `ciTools` to compute the two type of confidence intervals,
#' then we stack the dataframes together.
df1 <- add_ci(df, fit, names = c("lwr", "upr")) %>%
    mutate(type = "parametric")

df2 <- add_ci(df, fit, type = "boot", names = c("lwr", "upr"), nSims = 500) %>%
    mutate(type = "bootstrap")
df <- bind_rows(df1, df2)

#' This shows that there is a some difference between the parametric
#' method and the bootstrap method. However, the intervals disagree
#' the most where the estimated probability is close to $0$ or
#' $1$. When the estimated probability is moderate, we find that the
#' interval estimates essentially overlap.

#+ fig.width = 10, fig.heither = 7, fig.align = "center"
ggplot(df, aes(x = x, y = y)) +
    geom_jitter(height = 0.01) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
    geom_line(aes(x =x , y = pred), size = 2) +
    facet_grid(~type)

df3 <- filter(df, type == "parametric")

df4 <- filter(df, type == "bootstrap") %>%
    rename(lboot = lwr, uboot = upr) %>%
    bind_cols(df3)

#' Another perspective on the difference between these two interval
#' calculation methods. Notice that although the bootstrap intervals
#' are computed with an iterative method, the intervals contours are
#' smooth.

#+ fig.width = 9, fig.heither = 7, fig.align = "center"
ggplot(df4, aes(x = x, y = y)) +
    geom_jitter(height = 0.01) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill = "royalblue") +
    geom_ribbon(aes(ymin = lboot, ymax = uboot), alpha = 0.4, fill = "red") +
    geom_line(aes(x = x , y = pred), size = 2)
#'
#' If the sample increases, we will find that the two estimates
#' increasingly agree and converge to $0$ in width. Note that we do
#' not calculate prediction intervals for $y$ because the support of
#' $y$ is exactly $\{0,1\}$.
#'
#' ## Prediction Intervals
#'
#' Generally, parametric prediction intervals for GLMs are not
#' available. The solution `ciTools` takes is to perform a parametric
#' bootstrap on the model fit, then take quantiles on the bootstrapped
#' data produced for each observation. The procedure is performed
#' through the function `arm::sim`. The method of the parametric
#' bootstrap is described by the following algorithm:
#' 
#' 1. Fit the GLM, and collect the fitted responses. Set the number of
#' simulations, $M$.
#' 
#' 2. Simulate $M$ draws of the regression coefficients,
#' $\hat{\beta}_{*}$, from $N(\hat{\beta},
#' \hat{\mathrm{Cov}}(\hat{\beta}))$, where
#' $\hat{\mathrm{Cov}}(\hat{\beta}) = \hat{\sigma}^2 (X'X)^{-1}$.
#' 
#' 3. Simulate $[y_{*}|x]$ from the response distribution with mean
#' $g^{-1}(X \hat{\beta}_{*})$ and a variance determined by the
#' response distribution.
#' 
#' 4. Determine the $\alpha/2$ and $1-\alpha/2$ quantiles of the
#' simulated response $[y_{*}|x]$ for each $x$.
#'
#' This parametric bootstrap method propagates the uncertainty in the
#' regression effects $\hat{\beta}$ into the simulated draws from the
#' predictive distribution.
#'
#' Generally, there are many different ways to calculate the quantiles
#' of an empirical distribution, but the approach that `ciTools` takes
#' ensures that estimates quantiles lie in the support set of the
#' response. The choice we make corresponds to the `type = 1` argument
#' of `quantile()`.
#' 
#' We have seen in our simulations that the parametric bootstrap
#' provides (with some exception) interval estimates with
#' approximately nominal probability coverage. The unfortunate side
#' effect of opting to construct prediction intervals through a
#' parametric bootstrap is that the parameters of the predictive
#' distributions need to hard coded for each model. At this point, the
#' only models that we support with `add_pi()` are Guassian, Poisson,
#' Binomial, Gamma, and Quasipoisson (all of which need to be fit with
#' `glm()`).
#'
#' One exception to this scheme are GLMs with Gaussian errors. In this
#' case we may parametrically calculate prediction intervals. Under
#' Gaussian errors, `glm()` permits the use of the link functions
#' "identity", "log", and "inverse". The corresponding models are given
#' by the following:
#'
#' $$
#' \begin{equation}
#' \label{eq:gauss-link}
#' \begin{split}
#' y = X\beta + \epsilon\\
#' y = \exp(X\beta) + \epsilon \\
#' y = \frac{1}{X\beta} + \epsilon
#' \end{split}
#' \end{equation}
#' $$
#'
#' The conditional response distribution in each case may be written
#' $y \sim N(g^{-1}(X\beta), \sigma^2)$, and prediction intervals may be
#' computed parametrically. Parametric prediction intervals are therefore
#' constructed via
#'
#' $$
#' \begin{equation}
#' \label{eq:gauss-pi}
#' g^{-1}(x'\beta) \pm t_{1-\alpha/2, n-p-1} \sqrt{\hat{\sigma}^2 +
#'     \hat{\sigma}^2x'(X'X)^{-1}x} 
#' \end{equation}
#' $$
#'
#' Where, as in the linear model, $\hat{\sigma}^2$ estimates the
#' predictive uncertainty and $\hat{\sigma}^2 x(X'X)^{-1}x$ estimates
#' the uncertainty in the regression coefficients. Note that $g^{-1}(X
#' \hat{\beta})$ is a maximum likelihood estimate for the parameter
#' $g^{-1}(X \beta)$, by the functional invariance of maximum
#' likelihood estimation.
#' 
#' ## Poisson Example
#'
#' Poisson regression is usually the first line of defense against
#' count data, so we wish to present a complete example of quantifying
#' uncertainty for this type of model with `ciTools`. For simplicity
#' we fit a model on fake data.
#'
#'
#' We use rnorm to generate the a covariate, the randomness of $x$ has
#' no bearing on the model.
#' 
x <- rnorm(100, mean = 0)
y <- rpois(n = 100, lambda = exp(1.5 + 0.5*x))

df <- data.frame(x = x, y = y)
fit <- glm(y ~ x , family = poisson(link = "log"))

#' As seen previously, the commands in `ciTools` are "pipeable". Here,
#' we compute confidence and prediction intervals for a model fit at
#' the default level of $95\%$. The warning message only serves to
#' remind the user the precise quantiles cannot be formed for
#' non-continuous distributions.

df_ints <- df %>% 
    add_ci(fit, names = c("lcb", "ucb")) %>%
    add_pi(fit, names = c("lpb", "upb")) %>%
    print()

#' As with other methods available in `ciTools` the requested
#' statistics are computed, appended to the data frame, and returned to
#' the user as a tibble.

#+ fig.width = 10, fig.heither = 7, fig.align = "center"
df_ints %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point(size = 2) +
    geom_line(aes(x = x, y = pred), size = 1.2) + 
    geom_ribbon(aes(ymin = lcb, ymax = ucb), alpha = 0.4) +
    geom_ribbon(aes(ymin = lpb, ymax = upb), alpha = 0.2)

#' Since the response $y$ is count data, and the method we used to
#' determine the intervals involves simulation, we find that `ciTools`
#' will produce "jagged" bounds when all the intervals are plotted
#' simultaneously. 
#' 
#' We may also wish to compute response level probabilities and
#' quantiles. `ciTools` can also handle these estimates with
#' `add_probs()` and `add_quantile()` respectively. Since these
#' statistics are formed on the predictive distribution constructed by
#' `add_pi()`, we also estimates quantiles and probabilities using the
#' parametric bootstrap.
#'
df %>%
    add_probs(fit, q = 10) %>%
    add_quantile(fit, p = 0.4) %>%
    print()

#' ### Extension to Quasipoisson
#' 
#' A common problem with this model is the presence of
#' over-dispersion. Recall that for the Poisson model, we require that
#' the variance and mean agree, however this is practically a strict
#' and unreasonable modeling assumption. A quasipoisson model is one
#' remedy: it estimates a scale parameter as well and will provide a
#' better fit. Under quasipoisson assumption
#'
#' $$
#' E[y|X] = \mu = \exp (X \beta)
#' $$
#'
#' and
#'
#' $$
#' Var[y|X] = \phi \mu
#' $$
#'
#' Quasi models are not full maximum likelihood models, however it is possible
#' to embed a Quasipoisson in the Negative Binomial framework using
#'
#' $$
#' QP(\mu, \theta) = \mathrm{NB}(\mu, \theta = \frac{\mu}{\phi - 1})
#' $$
#'
#' Where NB is the parameterization of the Negative Binomial
#' distribution used by `glm.nb` in the `MASS` library. This model for
#' the negative binomial distribution, a continuous mixture of Poisson
#' random variables with gamma distributed means, is preferred over
#' that classical parameterization in applications. The preferences
#' stems from the fact that it allows for a real valued "$\theta$".
#'
#' Again, we generate fake data.
#' 
x <- runif(100, 0, 2)
mu <- exp(2 + 2 * x)
y <- rnegbin(n = 100, mu = mu, theta = mu/(6 - 1))

#' The data is over-dispersed: 
df <- data.frame(x = x, y = y)
fit <- glm(y ~ x, family = quasipoisson(link = "log"))
summary(fit)$dispersion

#' But `ciTools` can still construct appropriate intervals:
df_ints <- add_ci(df, fit, names = c("lcb", "ucb")) %>%
    add_pi(fit, names = c("lpb", "upb")) 

#+ fig.width = 10, fig.heither = 7, fig.align="center"
ggplot(df_ints, aes(x = x, y = y)) +
    geom_point(size = 2) +
    geom_line(aes(x = x, y = pred), size = 1.2) +
    geom_ribbon(aes(ymin = lcb, ymax = ucb), alpha = 0.4) +
    geom_ribbon(aes(ymin = lpb, ymax = upb), alpha = 0.2)


#' ## Other non-linear models
#'
#' Traditional negative binomial model fits are currently not
#' compatible with `ciTools`. While these models are similar to
#' Quasipoisson models, they require a different model fitting
#' algorithm than the IWLS method implemented in `glm`. Robust *t*
#' models can be fit with the `hett` package and are also not
#' supported by `ciTools` at this time. These models are fairly common
#' and we hope to support them in the future. 
