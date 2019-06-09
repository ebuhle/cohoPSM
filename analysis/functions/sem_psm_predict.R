#' Posterior Predictive Distribution of PSM Risk
#' 
#' Compute the posterior predictive distribution of PSM risk from a fitted SEM, 
#' possibly conditioned on specified values of the latent factor(s). 
#' 
#' Currently predictions can only be made for sites (with or without PSM data)
#' used to fit the model, and for sample-average precipitation conditions (i.e., precip effects
#' are set to zero). These restrictions will be lifted in a later version.
#'
#' @param fit Object of class \code{stanfit} representing a fitted PSM SEM. Parameters 
#' \code{mu_b0}, \code{sigma_b0}, \code{b0_std}, and \code{b0_Z} must be monitored.
#' @param data The data list passed to \code{stan()} to estimate \code{fit}.
#' @param newsites Numeric vector of length \code{N_new} giving the site numbers 
#' (coded as in \code{data}) for which predictions are to be made.
#' @param newZ Matrix of dimension \code{N_new x K} giving the \emph{K} factor scores
#' for each new prediction. If \code{NULL}, the posterior samples of site-specific factor
#' scores are used.
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

sem_psm_predict <- function(fit, data, newsites, newZ = NULL, level = "site", transform = FALSE) 
{
  samples <- extract(fit)
  
  logit_p_psm <- with(c(data, samples), {
    
    intercept <- as.vector(mu_b0) + as.vector(sigma_b0) * b0_std
    
    if(is.null(newZ)) {
      logit_p_psm_hat <- intercept[,newsites] + 
        I0_Z*apply(sweep(Z[,newsites,,drop = FALSE], c(1,3), b0_Z, "*"), c(1,2), sum) # add precip terms #
    } else {
      logit_p_psm_hat <- intercept[,newsites] + I0_Z*b0_Z  %*% t(newZ)  # add precip terms #
    }
    
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

