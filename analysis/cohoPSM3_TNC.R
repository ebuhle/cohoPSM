#==================================================================
# SETUP
#==================================================================

options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
library(yarrr)
library(vioplot)
library(colourvalues)
library(viridis)
library(shape)
library(rstan)
library(loo)
library(shinystan)
library(rstanarm)
# library(brms)
library(Hmisc)
library(matrixStats)
library(here)

# load functions
source(here("analysis","functions","stan_mean.R"))
source(here("analysis","functions","extract1.R"))
source(here("analysis","functions","stan_init.R"))
source(here("analysis","functions","stan_init_cv.R"))
source(here("analysis","functions","kfold_partition.R"))
source(here("analysis","functions","sem_psm_predict.R"))
source(here("analysis","functions","sem_z_crit.R"))
source(here("analysis","functions","vioplot2.R"))

# read and wrangle data
source(here("analysis","cohoPSM1_data.R"))  

# load previously saved stanfit objects
if(file.exists(here("analysis","results","glmm_psm.RData")))    # traffic GLMM
  load(here("analysis","results","glmm_psm.RData"))
if(file.exists(here("analysis","results","stan_psm.RData")))    # full SEM
  load(here("analysis","results","stan_psm.RData"))
if(file.exists(here("analysis","results","stan_psm_all.RData"))) # full SEM incl WADOE basins
  load(here("analysis","results","stan_psm_all.RData"))
if(file.exists(here("analysis","results","stan_psm_rt.RData"))) # roads + traffic SEM
  load(here("analysis","results","stan_psm_rt.RData"))
if(file.exists(here("analysis","results","stan_psm_tnc_cv_site.RData"))) # leave-site-out CV
  load(here("analysis","results","stan_psm_tnc_cv_site.RData"))


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
get_LL_glmm_psm

#==================================================================
# HIERARCHICAL REGRESSION MODELS FOR PSM 
#==================================================================

# Modify the PSM dataset to make it suitable for regression modeling.
# precip, traffic and log(traffic) centered and scaled to SD = 1 

## @knitr data_glmm_psm
psm_reg <- transform(psm, ID = 1:nrow(psm), ppt_su = scale(ppt_su), ppt_fa = scale(ppt_fa),
                     traffic = scale(traffic), log_traffic = scale(log(pmax(traffic, 0.1))))
## @knitr ignore

# Fit the "full" GLMM with summer and fall precip and log(traffic) as predictors.
# This model has the same structure as the "GLMM-like" submodel of the SEM, but 
# replaces the latent urbanization factor(s) with a single indicator variable, traffic.
# (The log transformation is used to emulate the log link function used for traffic and
# other gamma-distributed landscape variables in the "factor-analytic" submodel of the SEM.) 

## @knitr fit_glmm_psm
glmm_psm <- stan_glmer(cbind(n_psm, n - n_psm) ~ (ppt_su + ppt_fa) * log_traffic + 
                         (ppt_su + ppt_fa || site) + (1 | ID),
                       data = psm_reg,
                       family = binomial("logit"),
                       prior_intercept = normal(0,3),
                       prior = normal(0,3),
                       prior_covariance = decov(),
                       chains = 3, iter = 3000, warmup = 1000,
                       control = list(adapt_delta = 0.9))

# Inspect result
print(glmm_psm, digits = 2)
summary(glmm_psm, prob = c(0.025, 0.5, 0.975), pars = "beta", include = FALSE, digits = 2)

## @knitr ignore
# Save stanreg object
save(glmm_psm, file = here("analysis","results","glmm_psm.RData"))

# Calculate the marginal log-likelihood of the observed PSM frequencies under the GLMM.

## @knitr calc_LL_glmm_psm
LL_glmm_psm <- get_LL_glmm_psm(glmm_psm)
## @knitr ignore

#==================================================================
# STRUCTURAL EQUATION MODELS FOR LANDSCAPE ATTRIBUTES AND PSM
#==================================================================

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
stan_dat_rt <- list(S = nrow(X), 
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

## @knitr stan_psm_rt
# Fit it!
stan_psm_rt <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
                    data = stan_dat_rt, 
                    init = lapply(1:3, function(i) stan_init(stan_dat_rt)),
                    pars = c("a0","A","Z","phi","g_mu_X",
                             "mu_b0","b0_Z","sigma_b0","b0",
                             "mu_b_su","b_su_Z","sigma_b_su","b_su",
                             "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                             "sigma_psm","p_psm","ll_psm"), 
                    chains = 3, iter = 12000, warmup = 2000, thin = 5)

# Inspect and use shinystan to explore samples
print(stan_psm_rt, prob = c(0.025, 0.5, 0.975), 
      pars = c("b0","b_su","b_ppt_su","b_fa","b_ppt_fa","g_mu_X","p_psm","ll_psm","Z"), include = F)
# launch_shinystan(stan_psm_rt)

## @knitr ignore
# Save stanfit
save(stan_psm_rt, file = here("analysis","results","stan_psm_rt.RData"))

#==================================================================
# MODEL SELECTION
#==================================================================

#------------------------------------------------------------------
# In-sample model selection
#------------------------------------------------------------------

# Having fit all three candidate models to sites with PSM observations,
# now compare expected out-of-sample performance via PSIS-LOO

## @knitr loo_psm
loo_psm <- list(SEM_full = loo(stan_psm, pars = "ll_psm", 
                               r_eff = relative_eff(as.array(stan_psm, pars = "ll_psm"))),
                SEM_rt = loo(stan_psm_rt, pars = "ll_psm", 
                             r_eff = relative_eff(as.array(stan_psm_rt, pars = "ll_psm"))),
                GLMM = loo(LL_glmm_psm, 
                           r_eff = relative_eff(exp(LL_glmm_psm), 
                                                chain_id = rep(1:dim(as.array(glmm_psm))[2], 
                                                               each = dim(as.array(glmm_psm))[1]))))

compare(loo_psm$SEM_full, loo_psm$SEM_rt, loo_psm$GLMM)
compare(loo_psm$SEM_full, loo_psm$SEM_rt)
compare(loo_psm$SEM_full, loo_psm$GLMM)

#---------------------------------------------------------------------
# K-fold cross-validation over SITES:
# Leave out one or more sites of PSM data at a time,
# fit candidate models to training data and evaluate log posterior
# predictive density for the held-out observations.
# Sites are randomly partitioned into K = 10 groups that are
# roughly similar in size (i.e., number of observations).
#---------------------------------------------------------------------

# Randomly partition sites into groups
partitions <- kfold_partition(psm_dat = psm, K = 10)
partitions  # check that the procedure found a "good" partition (roughly equal group sizes)
site_group <- partitions$group

# List of candidate models
stan_psm_cv_site_list <- vector("list", 3)
names(stan_psm_cv_site_list) <- c("SEM_full","SEM_rt","GLMM")
for(i in 1:length(stan_psm_cv_site_list))
{
  stan_psm_cv_site_list[[i]] <- list(fit = vector("list",length(unique(site_group))),
                                     site_group = site_group,
                                     ll_psm = matrix(NA, 6000, nrow(psm)),
                                     elpd = NULL)
  names(stan_psm_cv_site_list[[i]]$fit) <- sort(unique(site_group))
}

# Loop over groups of sites, fit each candidate model to training data, 
# predict hold-out data, and store stanfit objects
for(j in sort(unique(site_group)))
{
  cat("Working on hold-out group ", grep(j, sort(unique(site_group))), "/", 
      length(unique(site_group)), " (see Viewer for progress) \n\n", sep = "")
  
  # Modify data to be passed to SEMs
  stan_dat_cv_site <- stan_dat
  stan_dat_cv_site$I_fit <- as.numeric(site_group != j)  # training data
  stan_dat_cv_site$I_lpd <- as.numeric(site_group == j)  # hold-out data
  stan_dat_rt_cv_site <- stan_dat_rt
  stan_dat_rt_cv_site$I_fit <- as.numeric(site_group != j)  # training data
  stan_dat_rt_cv_site$I_lpd <- as.numeric(site_group == j)  # hold-out data
  
  ## Full SEM
  # Fit
  fit <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
              data = stan_dat_cv_site, 
              init = lapply(1:3,function(i) stan_init_cv(stan_psm)),
              pars = c("a0","A","Z","phi","mu_b0","b0_Z","sigma_b0","b0",
                       "mu_b_su","b_su_Z","sigma_b_su","b_su",
                       "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                       "sigma_psm","p_psm","ll_psm"),
              chains = 3, iter = 12000, warmup = 2000, thin = 5,
              control = list(stepsize = 0.05))
  
  # Store fitted object and log-likelihood matrix
  stan_psm_cv_site_list$SEM_full$fit[[as.character(j)]] <- fit
  stan_psm_cv_site_list$SEM_full$ll_psm[,site_group == j] <- as.matrix(fit,"ll_psm")
  
  ## Roads + traffic SEM
  # Fit
  fit <- stan(file = here("analysis","stan","cohoPSM_SEM.stan"),
              data = stan_dat_rt_cv_site, 
              init = lapply(1:3,function(i) stan_init_cv(stan_psm_rt)),
              pars = c("a0","A","Z","phi","mu_b0","b0_Z","sigma_b0","b0",
                       "mu_b_su","b_su_Z","sigma_b_su","b_su",
                       "mu_b_fa","b_fa_Z","sigma_b_fa","b_fa",
                       "sigma_psm","p_psm","ll_psm"),
              chains = 3, iter = 12000, warmup = 2000, thin = 5,
              control = list(stepsize = 0.05))
  
  # Store fitted object and log-likelihood matrix
  stan_psm_cv_site_list$SEM_rt$fit[[as.character(j)]] <- fit
  stan_psm_cv_site_list$SEM_rt$ll_psm[,site_group == j] <- as.matrix(fit,"ll_psm")
  
  ## GLMM
  # Fit
  fit <- update(glmm_psm, subset = site_group != j)
  
  # Store fitted object and log-likelihood matrix
  stan_psm_cv_site_list$GLMM$fit[[as.character(j)]] <- fit
  ll <- get_LL_glmm_psm(fit, data = psm_reg[site_group == j,])
  stan_psm_cv_site_list$GLMM$ll_psm[,site_group == j] <- ll
}

# Average the posterior predictive density (not log density!) across posterior samples,
# take the log, then sum across observations to get expected log predictive density
for(i in 1:length(stan_psm_cv_site_list))
  stan_psm_cv_site_list[[i]]$elpd <- sum(log(colMeans(exp(stan_psm_cv_site_list[[i]]$ll_psm))))

# Summarize model comparison results
stan_psm_cv_site_mods <- data.frame(model = names(stan_psm_cv_site_list),
                                    elpd = sapply(stan_psm_cv_site_list,
                                                  function(x) x$elpd),
                                    se_elpd = sapply(stan_psm_cv_site_list, function(x) 
                                      sqrt(ncol(x$ll_psm))*sd(log(colMeans(exp(x$ll_psm))))),
                                    d_elpd = NA,
                                    se_d_elpd = NA)

for(i in 1:length(stan_psm_cv_site_list))
{
  stan_psm_cv_site_mods$d_elpd[i] <- stan_psm_cv_site_mods$elpd[i] - max(stan_psm_cv_site_mods$elpd)
  m <- which.max(stan_psm_cv_site_mods$elpd)
  devs <- log(colMeans(exp(stan_psm_cv_site_list[[i]]$ll_psm))) - 
    log(colMeans(exp(stan_psm_cv_site_list[[m]]$ll_psm)))
  stan_psm_cv_site_mods$se_d_elpd[i] <- sqrt(length(devs))*sd(devs)
  rm(m);rm(devs)
}

stan_psm_cv_site_mods

# Save objects
save(stan_psm_cv_site_list, stan_psm_cv_site_mods, file = here("analysis","results","stan_psm_tnc_cv_site.RData"))


#==================================================================
# FIGURES
#==================================================================

#----------------------------------------------------------------------
# PSM Threshold vs. Critical Urbanization Decision Analysis
#
# Plot site-specific posterior predictive distributions of PSM risk 
# as a function of Z, with posterior medians of current conditions.
# Show psm_crit and corresponding site-specific z_crit values for a
# specified confidence level alpha.
# Highlight a selected site and show detail
#----------------------------------------------------------------------

save_plot <- TRUE
psm_crit <- 0.3   # PSM threshold
alpha <- 0.9  # credibility level
prediction_level <- "site"  # "site" or "year"-within-site
show_site <- "400"
show_num <- which(levels(psm_all$site) == show_site)
show_crv <- ifelse(show_site %in% psm_all$site[psm_all$data=="psm"], show_num, stan_dat_all$S)

Z_draws <- extract1(stan_psm_all, "Z")[,,1]
Z <- colMedians(Z_draws)
PSM <- colMedians(sem_psm_predict(stan_psm_all, data = stan_dat_all, newsites = 1:stan_dat_all$S, 
                                  level = prediction_level, transform = TRUE))  # use estimated Z
z_out <- sem_z_crit(stan_psm_all, data = stan_dat_all, psm_crit = psm_crit, 
                    level = prediction_level, alpha = alpha)
newsites <- rep(sort(unique(c(stan_dat_all$site[stan_dat_all$I_fit==1], show_crv))), each = 50)
newZ <- matrix(rep(seq(min(Z, z_out$z_crit) - 0.1, 
                       max(Z, z_out$z_crit) + 0.1, 
                       length = 50), length(unique(newsites))), ncol = 1)
psm_pred <- sem_psm_predict(stan_psm_all, data = stan_dat_all, newsites = newsites, newZ = newZ,
                            level = prediction_level, transform = TRUE)
psm_pred_show_site <- sem_psm_predict(stan_psm_all, data = stan_dat_all, newsites = show_num,
                                      newZ = z_out$z_crit[show_num], transform = TRUE)

c1 <- transparent("darkgray", 0.3)  # posterior predictive PSM curves, all sites
c2 <- "gray"  # highlighted site
c2t <- transparent(c2, 0.6)
dzcols <- color_values(z_out$delta_z, palette = t(col2rgb(cividis(256, direction = -1))))
dzcolst <- transparent(dzcols, 0.1)

if(save_plot) {
  png(filename=here("analysis","results","figures","psm_z_threshold.png"),
      width=7.5, height=7, units="in", res=300, type="cairo-png") 
} else {
  dev.new(width = 7.5, height = 7)
}

par(mar = c(5.1, 4.2, 4.1, 4))
      
plot(Z, PSM, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlim = range(newZ), ylim = c(0,1), xaxs = "i",
     xlab = bquote("Urbanization (" * italic(Z) * ")"), ylab = "Predicted PSM")
# all PSM vs. Z curves and current conditions
for(j in unique(newsites))
  lines(newZ[newsites==j], colMedians(psm_pred[,newsites==j]), col = c1)
points(Z, PSM, pch = ifelse(psm_all$data=="psm", 16, 1), cex = 1.5, col = dzcolst)
# selected site: PSM vs. Z curve 
polygon(c(newZ[newsites==show_crv], rev(newZ[newsites==show_crv])),
        c(colQuantiles(psm_pred[,newsites==show_crv], probs = 0.05),
          rev(colQuantiles(psm_pred[,newsites==show_crv], probs = 0.95))),
        col = c2t, border = NA)
lines(newZ[newsites==show_crv], colMedians(psm_pred[,newsites==show_crv]), lwd = 3)
# selected site: posterior density of PSM at z_crit
vioplot2(psm_pred_show_site, at = z_out$z_crit[show_num], add = TRUE, 
         col = NULL, border = "black", lwd = 2, wex = 0.15, drawRect = FALSE, pchMed = "")
abline(v = z_out$z_crit[show_num], col = "red")
text(z_out$z_crit[show_num], par("usr")[3] - 0.03, bquote(italic(z)[crit]), adj = c(0.5,0), 
     xpd = TRUE, col = "red")
# selected site: delta_PSM (change in median PSM risk at z_crit relative to current)
segments(x0 = z_out$z_crit[show_num], x1 = Z[show_num], y0 = median(psm_pred_show_site), lty = 2)
arrows(x0 = Z[show_num], y0 = PSM[show_num], y1 = median(psm_pred_show_site), length = 0.1)
text(Z[show_num] + 0.02, 0.75*PSM[show_num] + 0.25*median(psm_pred_show_site), expression(Delta * "PSM"), adj = 0)
# selected site: posterior density of z; delta_z
vioplot2(Z_draws[,show_num], quantiles = alpha, horizontal = TRUE, at = PSM[show_num],
         add = TRUE, col = NULL, border = "black", lwd = 2, lwd.quantile = 2, wex = 0.05, 
         drawRect = FALSE, pchMed = "")
qz <- quantile(Z_draws[,show_num], alpha)
arrows(x0 = qz, x1 = z_out$z_crit[show_num], y0 = PSM[show_num], 
       col = dzcols[show_num], length = 0.1, lwd = 2)
text(0.25*qz + 0.75*z_out$z_crit[show_num], PSM[show_num] + 0.01, expression(Delta * italic(z)), 
     adj = c(0.5,0), col = dzcols[show_num])
# selected site: current conditions
points(Z[show_num], PSM[show_num], pch = ifelse(show_site %in% psm$site, 16, 1), 
       col = dzcols[show_num], cex = 1.5)
# PSM threshold and all z_crit values
abline(h = psm_crit, col = "red", lwd = 2)
text(par("usr")[1] - 0.02, psm_crit, bquote(PSM[crit]), adj = c(1,0.5), col = "red", xpd = TRUE)
rug(z_out$z_crit, col = transparent("red", 0.3))
shape::colorlegend(cividis(100, direction = -1, alpha = 0.9), 
                   zlim = round(range(z_out$delta_z)), dz = 1,
                   digit = 0, main = expression(Delta * italic(z)), main.cex = 1.5)

if(save_plot) dev.off()


#----------------------------------------------------------------------
# Delta z vs. Benefit
#
# Plot site-specific delta_z (calculated above) against a selected
# "benefit" (which may be delta_PSM or some independent attribute
# of subbasins). Divide plot into "restoration" (delta_z < 0)
# and "conservation" (delta_z > 0) regions and shade points by
# (standardized) benefit/delta_z, interpreted either as "return on
# investment" (for restoration) or "sensitivity" (for conservation).
#----------------------------------------------------------------------

benefit <- "delta_PSM"  # choose variable for y-axis

if(benefit == "delta_PSM") {
  psm_pred_z_crit <- sem_psm_predict(stan_psm_all, data = stan_dat_all, 
                                     newsites = 1:stan_dat_all$S, newZ = z_out$z_crit, 
                                     level = prediction_level, transform = TRUE)
  delta_psm <- colMedians(psm_pred_z_crit) - PSM
}

dbdzcols <- color_values(abs(delta_psm)/z_out$delta_z, palette = t(col2rgb(cividis(256, direction = -1))))
dbdzcolst <- transparent(dzcols, 0.1)
  
if(save_plot) {
  png(filename=here("analysis","results","figures","delta_z_vs_benefit.png"),
      width=7, height=7, units="in", res=300, type="cairo-png") 
} else {
  dev.new(width = 7, height = 7)
}

par(mar = c(5.1, 4.2, 4.1, 4))

plot(z_out$delta_z, delta_psm, pch = "", las = 1, cex.axis = 1.2, cex.lab = 1.5,
     xlab = expression(Delta * italic(z)), 
     ylab = ifelse(benefit=="delta_PSM", expression(Delta * "PSM"), benefit))
abline(v = 0)
points(z_out$delta_z, delta_psm, pch = 16, col = dbdzcolst, cex = 1.5)

if(save_plot) dev.off()



