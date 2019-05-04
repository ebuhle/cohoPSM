#+ setup_tnc
## @knitr setup_tnc
if(.Platform$OS.type == "windows") options(device=windows)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
library(rstan)
library(loo)
library(shinystan)
library(rstanarm)
library(brms)
library(Hmisc)
library(matrixStats)
library(here)
source(here("analysis","stan_mean.R"))
source(here("analysis","extract1.R"))
source(here("analysis","cohoPSM1_data.R"))  # read and wrangle data
# load previously saved stanfit objects and don't 
# re-evaluate the corresponding code chunks
if(file.exists(here("analysis","results","stan_psm.RData"))) {
  load(here("analysis","results","stan_psm.RData"))
  eval_stan_psm <- FALSE
}
if(file.exists(here("analysis","results","glmm_psm.RData"))) {
  load(here("analysis","results","glmm_psm.RData"))
  eval_glmm_psm <- FALSE
}
## @knitr ignore


#'==================================================================
#' HIERARCHICAL REGRESSION MODELS FOR LANDSCAPE DATA AND PSM
#'==================================================================

#+ data_psm_glmm
# Modify dataset
# precip, traffic and log(traffic) centered and scaled to SD = 1
## @knitr data_psm_glmm
psm_all_reg <- transform(psm_all, ppt_su = scale(ppt_su), ppt_fa = scale(ppt_fa),
                         traffic = scale(traffic), log_traffic = scale(log(pmax(traffic, 0.1))))

#+ fit_psm_glmm, eval = eval_psm_glmm
# Fit full model
## @knitr fit_psm_glmm
glmm_psm <- stan_glmer(cbind(n_psm, n - n_psm) ~ (ppt_su + ppt_fa) * log_traffic + 
                         (ppt_su + ppt_fa || site) + (1 | ID),
                       data = psm_all_reg, subset = data == "psm",
                       family = binomial("logit"),
                       prior_intercept = normal(0,3),
                       prior = normal(0,3),
                       prior_covariance = decov(),
                       chains = 3, iter = 2000, warmup = 1000,
                       control = list(adapt_delta = 0.9))

print(glmm_psm, digits = 2)
summary(glmm_psm, prob = c(0.025, 0.5, 0.975), pars = "beta", include = FALSE, digits = 2)

#' The following minimal example demonstrates how `posterior_linpred()`, and presumably
#' `posterior_predict()`, handle *new levels* of grouping factors that define group-varying
#' parameters. The documentation says that in this case the predictions "marginalize over the
#' relevant variables", but it's not clear whether that means marginalizing over the 
#' *hyperparameters* by drawing new $\MVN group-level parameters (as in `lme4::simulate` with 
#' `re.form = NA`), or marginalizing over the *group-level effects* themselves. Annoyingly
#' for us, it appears to be the latter. Let's have a look.
#+ posterior_linpred_reprex

dat <- data.frame(group = gl(10, 5), y = rnorm(10)[rep(1:10, each = 5)] + rnorm(50))
lmm <- stan_lmer(y ~ 1 + (1|group), data = dat, chains = 1, seed = 666)
lp1 <- posterior_linpred(lmm, newdata = data.frame(group = 9:12), re.form = NA)
head(lp1)
lp2 <- posterior_linpred(lmm, newdata = data.frame(group = 9:12), re.form = ~ (1|group))
head(lp2)
colSds(lp2[,3:4] - lp1[,3:4])
median(as.matrix(lmm, regex_pars = "Sigma"))

#' The first two columns in `newdata` correspond to groups present in the original sample,
#' while the second two correspond to previously unobserved groups. With `re.form = NA`, 
#' the group-level effects are set to `0` and the prediction for all groups is the hyper-mean 
#' (`(Intercept)`). It's not clear what `re.form = NULL` is doing. It's not using the hyper-mean,
#' but the two new groups are identical. Maybe it's drawing a *single* new group-level 
#' intercept at each posterior *sample*? Looks like I'll have to put this one to the Stan forums.
#' 
#' In the meantime, it's actually not that difficult to construct the new observation-level 
#' residuals by drawing from their hyperdistribution conditional on each sample of the hyper-mean
#' and hyper-SD. One approach would be to do this manually, adding the residuals
#' to the output of `posterior_linpred()` called with `re.form` leaving out the
#' observation-level term(s) in question. But that would entail figuring out the internal `lme4`
#' `Z`-matrix plumbing. A less efficient but much simpler brute-force approach is to just call
#' `posterior_predict()` repeatedly, once for each row in `newdata` with an unobserved group
#' (since a different draw from the hyperdistribution is returned each time).

#+ brms_fitted_newgroups
# Function to simulate from the posterior distribution of the linear predictor,
# including group-varying intercepts for new groups not included in the fitted data
## @knitr brms_fitted_newgroups
brms_fitted_newgroups <- function(object, newdata, re_formula = NULL) 
{
  lp <- matrix(NA, nrow(as.matrix(object)), nrow(newdata))
  re <- ranef(object)
  # rows of newdata that have new levels of any grouping factors
  newgroups <- sapply(names(re), function(f) !(newdata[,f] %in% rownames(re[[f]])))
  if(!all(newgroups))
    lp[,!newgroups] <- fitted(object, newdata = newdata[!newgroups,,drop = FALSE], 
                              re_formula = re_formula)
  if(any(newgroups))
    for(i in which(newgroups))
      lp[,i] <- fitted(object, newdata = newdata[i,,drop = FALSE], 
                             re_formula = re_formula, allow_new_levels = TRUE)
  return(lp)
}



#+ glmm_psm_loglik
# Calculate marginal log-likelihood, marginalizing over obs-level random errors
# (to do this, specify new levels of ID)
## @knitr glmm_psm_loglik
newdat <- subset(psm_all_reg, data == "psm")
newdat <- data.frame(idx = rep(1:nrow(newdat), each = 500), 
                     newdat[rep(1:nrow(newdat), each = 500),])
newdat$ID <- max(newdat$ID) + 1:nrow(newdat) # unobserved ID levels
### WTF? posterior_linpred() ignores ID-level intercept (but works with toy example)
lp <- posterior_linpred(glmm_psm, newdata = newdat)

ll_glmm_psm <- t(posterior_linpred(glmm_psm, newdata = newdat, transform = TRUE))
ll_glmm_psm <- apply(ll_glmm_psm, 1, function(p) dbinom(newdat$n_psm, newdat$n, p))
ll_glmm_psm <- aggregate(ll_glmm_psm, by = list(idx = newdat$idx), mean)


#' Let's take it from the top, this time using `brms` and `predict` instead of 
#' `rstanarm` and `posterior_linpred`
#+ brms_predict_reprex

dat <- data.frame(group = gl(10, 5), y = rnorm(10)[rep(1:10, each = 5)] + rnorm(50))
lmm <- brm(y ~ 1 + (1|group), data = dat, chains = 1)
lmm
lp1 <- fitted(lmm, newdata = data.frame(group = 9:12), re_formula = NA, summary = FALSE)
head(lp1)
lp2 <- fitted(lmm, newdata = data.frame(group = 9:12), re_formula = ~ (1|group), 
               allow_new_levels = TRUE, summary = FALSE)
head(lp2)
lp3 <- brms_fitted_newgroups(lmm, newdata = data.frame(group = 9:20), re_formula = ~ (1|group))
head(lp3)

