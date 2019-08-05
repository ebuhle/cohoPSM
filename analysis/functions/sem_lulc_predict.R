#' Posterior Predictive Distribution of Landscape Variables
#'
#' Compute the posterior predictive distribution of landscape variables from a
#' fitted SEM, possibly conditioned on specified values of the latent factor(s).
#'
#' @param fit Object of class \code{stanfit} representing a fitted PSM SEM.
#'   Parameters \code{a0}, \code{A}, and \code{phi} must be monitored.
#' @param data The data list passed to \code{stan()} to estimate \code{fit}.
#' @param d Integer in {1, ..., \emph{D}} indicating which landscape
#'   variable to simulate.
#' @param newZ Matrix of dimension \code{N_new x K} giving the \emph{K} factor
#'   scores for each new prediction. If \code{NULL}, the posterior samples of
#'   site-specific factor scores are used.
#' @param linpred Logical indicating whether to return the linear predictor (the
#'   default is \code{FALSE}).
#'
#' @return List with elements \describe{ \item{\code{X}} An \code{iter x
#'   N_new} matrix containing posterior samples of the predicted landscape variable.
#'   \item{\code{g_mu_X}} And \code{iter x N_new} matrix containing posterior samples of the linear
#'   predictor evaluated at the specified values of \emph{Z} }
#'
#' @export

sem_lulc_predict <- function(fit, data, d, newZ = NULL, linpred = FALSE)
{
  samples <- extract(fit)
  
  with(c(data, samples), {
    
    # linear predictor
    if(is.null(newZ)) {
      g_mu_X <- a0[,d] + apply(sweep(Z, c(1,3), A[,d,], "*"), 1:2, sum)
    } else {
      g_mu_X <- a0[,d] + tcrossprod(A[,d,], newZ)
    }
    
    # response distribution
    if(d <= D_normal) {  # normal
      X <- array(rnorm(prod(dim(g_mu_X)), g_mu_X, phi[,d]), dim(g_mu_X))
    } else {  # gamma
      X <- array(rgamma(dim(g_mu_X), phi[,d], phi[,d] * exp(-g_mu_X)), dim(g_mu_X))
    }

    list(g_mu_X = if(linpred) g_mu_X else NULL, X = X)
  })
}
