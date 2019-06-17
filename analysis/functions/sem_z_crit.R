#' Critical Z-Values from an SEM Given a PSM Threshold
#'
#' Find site-specific critical values of a single latent factor \emph{Z} in a fitted SEM 
#' corresponding to a given critical PSM risk, and differences between the critical values
#' and current values of \emph{Z}, for a specified confidence level.
#' 
#' \code{sem_z_crit} assumes \code{fit} represents a single-factor model in which 
#' \code{I_b0_Z == 1} (that is, the site-specific fitted intercept includes an effect of \emph{Z},
#' and that the estimate of that effect is generally positive. Currently 
#' \code{sem_z_crit} assumes sample-average precipitation conditions; this will be relaxed in
#' the future. Critical values of \emph{Z} are found by first solving for the site-specific
#' zero of \emph{P(PSM | Z, data) = 0} for each posterior sample, and then finding the 
#' \emph{(1 - alpha/(1 - F(0))} quantile of those sampled zeros for each site that correspond to
#' positive slopes, where \emph{F(0)} is the empirical CDF of the slope w.r.t. \emph{Z} 
#' evaluated at zero. This method is conservative; it guarantees that at least proportion 
#' \emph{alpha} of the posterior sample of curves (i.e., those with positive slopes) will be 
#' less than \code{psm_crit} for all \emph{Z < Z_crit}, but it requires that 
#' \emph{F(0) < 1 - alpha}. Otherwise an error will be returned. Differences between \emph{z_crit} 
#' and the \emph{alpha}-th posterior quantile of estimated site-specific \emph{Z} are also
#' calculated.
#'
#' @param fit Object of class \code{stanfit} representing a fitted PSM SEM. Parameters 
#' \code{mu_b0}, \code{sigma_b0}, \code{b0_std}, and \code{b0_Z} must be monitored.
#' @param data The data list passed to \code{stan()} to estimate \code{fit}.
#' @param psm_crit PSM threshold in (0,1).
#' @param level Level of grouping at which to predict. Options are \code{"site"}
#'  (site-specific interannual average, the default) or \code{"year"}
#'  (include year-within-site residual variation).
#' @param alpha Confidence level.
#' 
#' @return A list with elements \describe{ 
#' \item{\code{zeros}} An iter x \emph{S} matrix, 
#' where \emph{S} is the number of sites in \code{fit}, containing posterior samples of 
#' the site-specific zeros of \emph{P(PSM | Z) - psm_crit} as a function of \emph{Z}
#' \item{\code{z_crit}} Numeric vector of length \emph{S} containing site-specific critical
#' \emph{Z}-values.
#' \item{\code{delta_z}} Numeric vector of length \emph{S} containing differences between
#' \code{z_crit} and the \emph{alpha}-th posterior quantile of site-specific \emph{Z}.
#' }
#' 
#' @export

sem_z_crit <- function(fit, data, psm_crit, level = "site", alpha) 
{
  samples <- extract(fit)
  logit_psm_crit <- qlogis(psm_crit)
  
  out <- with(c(samples, data), {
    
    iter <- nrow(b0_std)
    S <- ncol(b0_std)
    
    intercept <- as.vector(mu_b0) + as.vector(sigma_b0) * b0_std # add precip effects #
    if(level == "year") {
      delta <- matrix(rnorm(iter*S, 0, sigma_psm), iter, S)
      intercept <- intercept + delta
    }
    
    slope <- b0_Z  # add precip x Z interactions #
    cdf <- ecdf(slope)
    F0 <- cdf(0)
    if(F0 >= 1 - alpha) 
      stop(paste0("alpha = ", alpha, ", but P(slope > 0 | data) = ",  F0))
    
    zeros <- sweep(logit_psm_crit - intercept, 1, slope, "/")
    z_crit <- matrixStats::colQuantiles(zeros[slope > 0,], probs = 1 - alpha / (1 - F0))
    delta_z <- z_crit - matrixStats::colQuantiles(Z[,,1], probs = alpha)

    list(zeros = zeros, z_crit = z_crit, delta_z = delta_z)
  })
  
  return(out)
}