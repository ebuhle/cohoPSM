#------------------------------------------------------
# Function to generate initial values for SEM
#------------------------------------------------------

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

