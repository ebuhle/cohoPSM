#==================================================================
# SETUP
#==================================================================

if(.Platform$OS.type == "windows") options(device=windows)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
library(yarrr)
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
load(here("analysis","stan_init.RData"))
load(here("analysis","stan_init_cv.RData"))
load(here("analysis","KfoldCV_partition.RData"))
source(here("analysis","cohoPSM1_data.R"))  # read and wrangle data
# load previously saved stanfit objects
if(file.exists(here("analysis","results","stan_psm.RData"))) # full SEM from Roads to Ruin paper
  load(here("analysis","results","stan_psm.RData"))
if(file.exists(here("analysis","results","glmm_psm.RData")))
  load(here("analysis","results","glmm_psm.RData"))


#=================================================================================
# CALCULATE MARGINAL LIKELIHOOD FOR A GLMM BY INTEGRATING OUT RANDOM EFFECTS
#=================================================================================

#--------------
# rstanarm
#--------------

# The following minimal example demonstrates how rstanarm::posterior_linpred(), and presumably
# rstanarm::posterior_predict(), handle *new levels* of grouping factors that define group-varying
# parameters. The documentation says that in this case the predictions "marginalize over the
# relevant variables", but it's not clear whether that means marginalizing over the 
# *hyperparameters* by drawing new MVN group-level parameters (as in lme4::simulate with 
# re.form = NA), or marginalizing over the *group-level effects* themselves. Annoyingly
# for us, it appears to be the latter. Let's have a look.

set.seed(123)
dat <- data.frame(group = gl(10, 5), y = rnorm(10)[rep(1:10, each = 5)] + rnorm(50))
lmm <- stan_lmer(y ~ 1 + (1|group), data = dat, chains = 1, seed = 456)
lmm
lp_fixef <- posterior_linpred(lmm, newdata = data.frame(group = 9:12), re.form = NA)
head(lp_fixef)
lp_ranef <- posterior_linpred(lmm, newdata = data.frame(group = 9:12), re.form = ~ (1|group))
head(lp_ranef)
lp_ranef <- posterior_linpred(lmm, newdata = data.frame(group = 9:12), re.form = ~ (1|group))
head(lp_ranef)

# The first two columns in newdata correspond to groups present in the original sample,
# while the next two correspond to previously unobserved groups. With re.form = NA, 
# the group-level effects are set to 0 and the prediction for all groups is the hyper-mean 
# ((Intercept)). It's not clear what re.form = ~ (1|group) is doing. It's not using the hyper-mean,
# but the two new groups are identical, and repeated calls give the same result. Maybe it's 
# drawing a *single* new group-level intercept at each posterior *sample*? Looks like I'll 
# have to put this one to the Stan forums.

#--------------
# brms
#--------------

# Let's take it from the top, this time using brms and predict() instead of 
# rstanarm and posterior_linpred()

set.seed(123)
dat <- data.frame(group = gl(10, 5), y = rnorm(10)[rep(1:10, each = 5)] + rnorm(50))
lmm <- brm(y ~ 1 + (1|group), data = dat, chains = 1, seed = 456)
lmm
lp_fixef <- fitted(lmm, newdata = data.frame(group = 9:12), re_formula = NA, summary = FALSE)
head(lp_fixef)
lp_ranef <- fitted(lmm, newdata = data.frame(group = 9:12), re_formula = ~ (1|group), 
                   allow_new_levels = TRUE, summary = FALSE)
head(lp_ranef)
lp_ranef <- fitted(lmm, newdata = data.frame(group = 9:12), re_formula = ~ (1|group), 
                   allow_new_levels = TRUE, summary = FALSE)
head(lp_ranef)

# OK, so apparently predict() is generating new groups' values from the hyperdistribution (good)
# but reusing the same values for every new group (bad). However, it does appear to draw new values
# each time it's called, which suggests we could hack it by calling it repeatedly over all 
# the new groups. Except...there's something fishy about that second call. The first new value is
# suspiciously similar to that of group 10 (column 2). Let's try a few more, focusing on just a
# single new group.

lp_ranef <- sapply(1:5, function(i) fitted(lmm, newdata = data.frame(group = 11), re_formula = ~ (1|group), 
                                           allow_new_levels = TRUE, summary = FALSE))
head(lp_ranef)

# What's going on here? Several of these values are almost identical to each other and/or to
# values returned above for previously observed groups. Am I seeing things? Nope,
# a pairs plot of these "randomly generated" values shows some distinctly nonrandom pattern:
 
pairs(lp_ranef, col = transparent("darkblue", 0.7), labels = paste("call", 1:5))

#------------------------
# Kludgy Solutions
#------------------------

# Unless and until this is resolved, it looks like we're stuck with generating the linear predictor
# for new observation-level random effects by brute force. We can do this in rstanarm by adding 
# residuals drawn from the hyperdistribution to the output of posterior_linpred() called with a 
# re.form that leaves out the observation-level term(s) in question. This kludgy solution is limited
# to the case of random effects on the intercept (vs. the slopes).

# Function to simulate from the posterior distribution of the linear predictor,
# including group-varying intercepts for new groups not included in the fitted data
posterior_linpred_newgroups <- function(object, newdata) 
{
  lp <- matrix(NA, nrow(as.matrix(object)), nrow(newdata))
  re <- ranef(object)
  re_form <- lme4::findbars(formula(object))  # all random effect terms
  # rows of newdata that have new levels of any grouping factors
  newgroups <- sapply(names(re), function(f) {
    ng <- !(newdata[,f] %in% rownames(re[[f]]))
    if(any(ng) & any(colnames(re[[f]]) != "(Intercept)"))
      stop("New groups allowed only allowed for term (Intercept)") else {
        ## won't work with sapply ## re_form <- re_form - formula(paste0("~ (1|", f, ")"))
        return(ng)
      }
    })
  if(!all(newgroups))
    lp[,!newgroups] <- posterior_linpred(object, newdata = newdata[!newgroups,,drop = FALSE])
  if(any(newgroups)) {
    lp[,newgroups] <- posterior_linpred(object, newdata = newdata[i,,drop = FALSE], 
                                        re.form = re_form)
    ## ...and then somehow add random effects generated from the hyperdistributions corresponding
    ## to the discarded group-level intercepts???
  }
  return(lp)
}

# Arrgh, this is too convoluted for what is still a hacky and inflexible solution. 
# Might as well just tailor it specifically for the PSM GLMM case.
# The following function calculates the marginal log-likelihood of observed PSM frequencies
# from a fitted stanreg object by Monte Carlo integration over the observation-level residuals.

## @knitr get_LL_glmm_psm
get_LL_glmm_psm <- function(object, data = NULL, N_MC = 1000)
{
  if(is.null(data)) data <- object$data
  # random effects-only formula
  re_formula <- reformulate(sapply(lme4::findbars(formula(object)),
                                   function(f) paste("(", deparse(f), ")")))
  # drop obs-level random effect (assumes re_formula contains other random terms)
  re_formula <- update(re_formula, ~ . - (1|ID))
  # posterior draws of linear predictor without obs-level random effect
  lp <- posterior_linpred(object, newdata = data, re.form = re_formula)
  # generate new obs-level residuals, calculate log of marginal likelihood
  sigma_psm <- as.matrix(object, regex_pars= "Sigma\\[ID")
  LL <- sapply(1:nrow(data), function(i) {
    resid_psm_mc <- matrix(rnorm(nrow(sigma_psm)*N_MC, 0, sigma_psm), nrow(sigma_psm), N_MC)
    p_psm_mc <- plogis(lp[,i] + resid_psm_mc)
    LL_psm_mc <- dbinom(data$n_psm[i], data$n[i], p_psm_mc, log = TRUE)
    return(matrixStats::rowLogSumExps(LL_psm_mc) - log(N_MC))
  })
}
## @knitr ignore

#==================================================================
# HIERARCHICAL REGRESSION MODELS FOR PSM 
#==================================================================

# Modify the PSM dataset to make it suitable for regression modeling.
# precip, traffic and log(traffic) centered and scaled to SD = 1 

## @knitr data_glmm_psm
psm_all_reg <- transform(psm_all, ppt_su = scale(ppt_su), ppt_fa = scale(ppt_fa),
                         traffic = scale(traffic), log_traffic = scale(log(pmax(traffic, 0.1))))
## @knitr ignore

# Fit the "full" GLMM with summer and fall precip and log(traffic) as predictors.
# This model has the same structure as the "GLMM-like" submodel of the SEM, but 
# replaces the latent urbanization factor(s) with a single indicator variable, `traffic`.
# (The log transformation is used to emulate the log link function used for `traffic` and
# other gamma-distributed landscape variables in the "factor-analytic" submodel of the SEM.) 

## @knitr fit_glmm_psm
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
## @knitr ignore

# Calculate the marginal log-likelihood of the observed PSM frequencies under the GLMM.

## @knitr calc_LL_glmm_psm
LL_glmm_psm <- get_LL_glmm_psm(glmm_psm, data = psm_all_reg[psm_all_reg$data=="psm",])
## @knitr ignore

#==================================================================
# STRUCTURAL EQUATION MODELS FOR LANDSCAPE ATTRIBUTES AND PSM
#==================================================================

#------------------------------------------------------
# Function to generate initial values for SEM
#------------------------------------------------------

## @knitr stan_init
# Function to generate initial values for chains
stan_init <- function(stan_dat) 
{
  with(stan_dat, {
    # X <- t(X)
    D <- D_normal + D_gamma
    
    list(a0 = rnorm(D, c(colMeans(X[,0:D_normal]), colMeans(log(X[,(D_normal+1):D]))), 1),
         A_nid_vec = array(rnorm(D*L - L*(L-1)/2, 0, 1), dim = D*L - L*(L-1)/2),
         Z_nid = matrix(rnorm(S*L, 0, 1), nrow = S, ncol = L),
         phi = runif(D, 0.5, 1),
         # mu_b0 = rnorm(1, qlogis(mean(n_psm/n, na.rm=T)), 1),
         mu_b0 = rnorm(1, -3, 1),
         b0_Z_nid = array(rnorm(L, 1, 0.5), dim = L),
         sigma_b0 = runif(1, 1, 2),
         b0_std = array(rnorm(S, 0, 0.1), dim = S),
         mu_b_su = rnorm(1, 0.3, 0.1),
         b_su_Z_nid = array(rnorm(L, -0.2, 0.1), dim = L),
         sigma_b_su = runif(1, 0.1, 1),
         b_su_std = array(rnorm(S, 0, 0.1), dim = S),
         mu_b_fa = rnorm(1, 0.03, 0.02),
         b_fa_Z_nid = array(rnorm(L, -0.03, 0.02), dim = L),
         sigma_b_fa = runif(1, 0.02, 0.1),
         b_fa_std = array(rnorm(S, 0, 0.1), dim = S),
         sigma_psm = runif(1, 0.5, 2),
         logit_p_psm_std = array(rnorm(N, 0, 0.1), dim = N))
  })
}
## @knitr ignore

#------------------------------------------------------
# Model with roads and traffic only
#------------------------------------------------------

#------------------------------------------------------
# Assemble data for sites with PSM observations
# in Stan-friendly format
#
# Note that level 3 roads are omitted because
# including them leads to nontrivial multimodality
# and horrible mixing 
#------------------------------------------------------

## @knitr stan_data_psm_roads
# nonnegative continuous data:
# scale to SD = 1, bound away from 0
X <- as.matrix(lulc_roads_data[,c("roads1","roads2","roads4","roads5","traffic")])
X <- sweep(X, 2, apply(X, 2, sd), "/")
X[X==0] <- 1e-4
X <- sweep(X, 2, apply(X, 2, sd), "/")

# all covariates assumed to be gamma distributed
normal_indx <- NULL
gamma_indx <- 1:ncol(X)

# Data for Stan
stan_dat_roads <- list(S = nrow(X), 
                       D_normal = length(normal_indx), D_gamma = length(gamma_indx),
                       X = X, 
                       L = 1,  # user-specified!
                       N = nrow(psm), 
                       site = as.numeric(psm$site),
                       ppt_su = array(as.vector(scale(psm$ppt_su/10, scale = FALSE)), dim = nrow(psm)),
                       ppt_fa = array(as.vector(scale(psm$ppt_fa/10, scale = FALSE)), dim = nrow(psm)),
                       I0_Z = 1,
                       I_su = 1,
                       I_su_Z = 1,
                       I_fa = 1,
                       I_fa_Z = 1,
                       n = psm$n,
                       n_psm = psm$n_psm,
                       I_fit = rep(1, nrow(psm)),
                       I_lpd = rep(1, nrow(psm)))
## @knitr ignore

#-------------------------------------------------------------------
# Fit SEM with roads and traffic to sites with PSM observations
#-------------------------------------------------------------------

## @knitr stan_psm_roads
# Fit it!
stan_psm_roads <- stan(file = here("analysis","cohoPSM_SEM.stan"),
                       data = stan_dat_roads, 
                       init = lapply(1:3, function(i) stan_init(stan_dat_roads)),
                       pars = c("a0","A","Z","phi","g_mu_X",
                                "mu_b0","b0_Z","sigma_b0","b0",
                                "mu_b_su","b_su_Z","sigma_b_su","b_su",
                                "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                                "sigma_psm","p_psm","ll_psm"), 
                       chains = 3, iter = 12000, warmup = 2000, thin = 5)

# Inspect and use shinystan to explore samples
print(stan_psm_roads, prob = c(0.025, 0.5, 0.975), 
      pars = c("b0","b_su","b_ppt_su","b_fa","b_ppt_fa","g_mu_X","p_psm","ll_psm","Z"), include = F)
# launch_shinystan(stan_psm_roads)

# Save stanfit
save(stan_psm_roads, file = here("analysis","results","stan_psm_roads.RData"))
## @knitr ignore

#==================================================================
# Model Selection
#==================================================================

#------------------------------------------------------------------
# In-sample model selection
#------------------------------------------------------------------

# Having fit all three candidate models to sites with PSM observations,
# now compare expected out-of-sample performance via PSIS-LOO

## @knitr loo_psm
loo_psm <- list(stan_psm = loo(stan_psm, pars = "ll_psm", 
                               r_eff = relative_eff(as.array(stan_psm, pars = "ll_psm"))),
                stan_psm_roads = loo(stan_psm_roads, pars = "ll_psm", 
                               r_eff = relative_eff(as.array(stan_psm_roads, pars = "ll_psm"))),
                glmm_psm = loo(LL_glmm_psm, 
                               r_eff = relative_eff(exp(LL_glmm_psm), 
                                                    chain_id = rep(1:dim(as.array(glmm_psm))[2], 
                                                                   each = dim(as.array(glmm_psm))[1]))))

compare(loo_psm$stan_psm, loo_psm$stan_psm_roads)
compare(loo_psm$stan_psm, loo_psm$glmm_psm)

# ## K-fold cross-validation over SITES

# Leave out one or more sites of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations.
# Sites are randomly partitioned into K = 10 groups that are
# roughly similar in size (i.e., number of observations).









