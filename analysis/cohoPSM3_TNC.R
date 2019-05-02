## @knitr setup_tnc
if(.Platform$OS.type == "windows") options(device=windows)
library(rstan)
library(loo)
library(shinystan)
library(rstanarm)
library(Hmisc)
library(here)
source(here("analysis","stan_mean.R"))
source(here("analysis","extract1.R"))
source(here("analysis","cohoPSM1_data.R"))  # read and wrangle data
# load previously saved stanfit object containing "full" SEM
# and don't re-evaluate that code chunk
if(file.exists(here("analysis","results","stan_psm.RData"))) {
  load(here("analysis","results","stan_psm.RData"))
  eval_stan_psm <- FALSE
}
## @knitr ignore


#==================================================================
# HIERARCHICAL REGRESSION MODELS FOR LANDSCAPE DATA AND PSM
#==================================================================

# Modify dataset
# precip, traffic and log(traffic) centered and scaled to SD = 1
## @knitr data_psm_glmm
psm_all_reg <- transform(psm_all, ppt_su = scale(ppt_su), ppt_fa = scale(ppt_fa),
                         traffic = scale(traffic), log_traffic = scale(log(pmax(traffic, 0.1))))

# Fit full model
## @knitr fit_psm_glmm
glmm_psm <- stan_glmer(cbind(n_psm, n - n_psm) ~ (ppt_su + ppt_fa) * log_traffic + 
                         (ppt_su + ppt_fa || site) + (1 | ID),
                       data = psm_all_reg, subset = data == "psm",
                       family = binomial("logit"),
                       prior_intercept = normal(0,3),
                       prior = normal(0,3),
                       prior_covariance = decov(),
                       chains = 3, cores = 3, iter = 2000, warmup = 1000,
                       control = list(adapt_delta = 0.9))

print(glmm_psm, digits = 2)
summary(glmm_psm, prob = c(0.025, 0.5, 0.975), pars = "beta", include = FALSE, digits = 2)

# Calculate marginal log-likelihood, marginalizing over obs-level random errors
# (to do this, specify new levels of ID)
newdat <- psm_all_reg[psm_all_reg$data == "psm",]
newdat <- data.frame(idx = rep(1:nrow(newdat), each = 500), 
                     newdat[rep(1:nrow(newdat), each = 500),])
newdat$ID <- max(newdat$ID) + 1:nrow(newdat)
### WTF? posterior_linpred() ignores ID-level intercept (but works with toy example)
ll_glmm_psm <- t(posterior_linpred(glmm_psm, newdata = newdat, transform = TRUE))
ll_glmm_psm <- apply(ll_glmm_psm, 1, function(p) dbinom(newdat$n_psm, newdat$n, p))
ll_glmm_psm <- aggregate(ll_glmm_psm, by = list(idx = newdat$idx), mean)

### just for kicks
load(here("results","stan_psm.RData"))  # saved stanfit object from next code chunk below
psm_all_reg$Z[psm_all_reg$data == "psm"] <- stan_mean(stan_psm,"Z")[as.numeric(psm$site)]

glmm_psm_Z <- stan_glmer(cbind(n_psm, n - n_psm) ~ (ppt_su + ppt_fa) * Z + 
                           (ppt_su + ppt_fa || site) + (1 | ID),
                         data = psm_all_reg, subset = data == "psm",
                         family = binomial("logit"),
                         prior_intercept = normal(0,3),
                         prior = normal(0,3),
                         prior_covariance = decov(),
                         chains = 3, cores = 3, iter = 2000, warmup = 1000,
                         control = list(adapt_delta = 0.9))

glmm_psm_Z
summary(glmm_psm_Z, prob = c(0.025, 0.5, 0.975), pars = "beta", include = FALSE)
###


