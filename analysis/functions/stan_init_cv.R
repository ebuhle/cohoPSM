#---------------------------------------------------------------------
# Function to generate initial values for SEM for cross-validation,
# using previously fitted full model 
#---------------------------------------------------------------------

stan_init_cv <- function(fit)
{
  samples <- extract(fit)
  M <- nrow(samples$a0)
  i <- sample(M,1)
  
  AA <- matrix(samples$A[i,,], nrow = dim(samples$A)[2], ncol = dim(samples$A)[3], byrow = TRUE)
  AA <- AA[lower.tri(AA, diag = TRUE)]
  logit_p_psm_hat <- stan_mean(fit,"logit_p_psm_hat")
  logit_p_psm <- qlogis(samples$p_psm[i,])
  
  with(samples, 
       list(a0 = a0[i,],
            A_nid_vec = array(AA, dim = length(AA)),
            Z_nid = matrix(Z[i,,], nrow = dim(Z)[2], ncol = dim(Z)[3], byrow = TRUE),
            phi = phi[i,],
            mu_b0 = mu_b0[i],
            b0_Z_nid = array(b0_Z[i,], dim = ncol(b0_Z)),
            sigma_b0 = sigma_b0[i],
            b0_std = array((b0[i,] - mu_b0[i]) / sigma_b0[i], dim = ncol(b0)),
            mu_b_su = mu_b_su[i],
            b_su_Z_nid = array(b_su_Z[i,], dim = ncol(b_su_Z)),
            sigma_b_su = sigma_b_su[i],
            b_su_std = array((b_su[i,] - mu_b_su[i]) / sigma_b_su[i], dim = ncol(b_su)),
            mu_b_fa = mu_b_fa[i],
            b_fa_Z = array(b_fa_Z[i,], dim = ncol(b_fa_Z)),
            sigma_b_fa = sigma_b_fa[i],
            b_fa_std = array((b_fa[i,] - mu_b_fa[i]) / sigma_b_fa[i], dim = ncol(b_fa)),
            sigma_psm = sigma_psm[i],
            logit_p_psm_std = array((logit_p_psm - logit_p_psm_hat) / sigma_psm[i], dim = length(logit_p_psm))))
}
