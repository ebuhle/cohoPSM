# Compute predicted values of PSM risk from a fitted SEM, conditioned on input values of the
# latent factor(s). Currently predictions can only be made for sites (with or without PSM data)
# used to fit the model, and for sample-average precipitation conditions (i.e., precip effects
# set to zero). These restrictions will be lifted in a later version.
#
# fit   Object of class `stanfit` representing a fitted PSM SEM. arameters `b0`, `b0_Z`,
#       `b_ppt_su`, `b_su_Z`, `b_ppt_fa`, `b_fa_Z` must be monitored.
# data  The data list passed to `stan()` to estimate `fit`.
# newsites  A numeric vector of length `N_new` giving the site numbers (coded as in `data`) for which
#           predictions are to be made.
# newZ  A matrix of dimension N_new x K giving the K factor scores for each new prediction.
#
## Value
# An `iter` x `N_new` matrix containing posterior samples of  the predicted probability 
# (or logit probability, if `transform = FALSE`) of PSM for each 
# 

sem_psm_predict <- function(fit, data, newsites, newZ) 
{
  samples <- extract(fit)
  
  
}

