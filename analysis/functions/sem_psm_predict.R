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
# transform  Logical indicating whether to return the linear predictor (FALSE, the default) or
#            inverse-logit transform it.
#
## Value
# An `iter` x `N_new` matrix containing posterior samples of the predicted probability 
# (or the linear predictor of logit probability, if `transform = FALSE`) of PSM for each 
# 

sem_psm_predict <- function(fit, data, newsites, newZ, transform = FALSE) 
{
  samples <- extract(fit)

  with(c(data, samples), {
    logit_p_psm_hat <- b0 + I0_Z*b0_Z  %*% t(newZ) + 
      I_su*b_su + I_su_Z*b_su_Z %*% t(newZ) + I_fa_Z*b_fa_Z %*% t(newZ)
  })
  
  return(switch(transform, "FALSE" = logit_p_psm_hat, "TRUE" = plogis(logit_p_psm_hat)))
}

