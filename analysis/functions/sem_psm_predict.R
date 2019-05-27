#' Posterior Predictive Distribution of PSM Risk Given Factor Scores
#' 
#' Compute the posterior predictive distribution of PSM risk from a fitted SEM, 
#' conditioned on input values of the latent factor(s). 
#' 
#' Currently predictions can only be made for sites (with or without PSM data)
#' used to fit the model, and for sample-average precipitation conditions (i.e., precip effects
#' are set to zero). These restrictions will be lifted in a later version.
#'
#' @param fit Object of class \code{stanfit} representing a fitted PSM SEM. Parameters 
#' \code{b0}, \code{b0_Z}, \code{b_ppt_su}, \code{b_su_Z}, \code{b_ppt_fa}, \code{b_fa_Z}
#' must be monitored.
#' @param data The data list passed to \code{stan()} to estimate \code{fit}.
#' @param newsites Numeric vector of length \code{N_new} giving the site numbers 
#' (coded as in \code{data}) for which predictions are to be made.
#' @param newZ Matrix of dimension \code{N_new x K} giving the \emph{K} factor scores
#' for each new prediction.
#' @param level Level of grouping at which to predict. Options are \code{"site"}
#' (site-specific interannual average, the default) or \code{"year"} 
#' (include year-within-site residual variation).
#' @param transform Logical indicating whether to return the linear predictor (\code{FALSE}, 
#' the default) or inverse-logit transform it.
#'
#' @return An \code{iter x N_new} matrix containing posterior samples of the predicted probability 
#' (or the linear predictor of logit probability, if \code{transform == FALSE}) of PSM. 
#' 
#' @export

sem_psm_predict <- function(fit, data, newsites, newZ, level = "site", transform = FALSE) 
{
  samples <- extract(fit)
  
  logit_p_psm <- with(c(data, samples), {
    logit_p_psm_hat <- b0 + I0_Z*b0_Z  %*% t(newZ)  # add precip terms
    if(level == "site")
      logit_p_psm <- logit_p_psm_hat
    if(level == "year") {
      delta <- array(rnorm(prod(dim(logit_p_psm_hat)), 0, sigma_psm), dim(logit_p_psm_hat))
      logit_p_psm <- logit_p_psm_hat + delta
    }
    logit_p_psm
  })

  return(if(transform) plogis(logit_p_psm) else logit_p_psm)
}

