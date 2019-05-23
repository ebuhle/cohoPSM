data {
  int<lower=1> S;                // number of sites
  int<lower=0> D_normal;         // dimension of normally distributed site-level covariates
  int<lower=0> D_gamma;          // dimension of gamma-distributed site-level covariates
  matrix[S,D_normal+D_gamma] X;  // site-level covariates: 1, ..., D_normal, ..., D_normal+D_gamma
  int<lower=1> L;                // dimension of latent factor space (user-specified!)
  int<lower=1> N;                // number of PSM observations
  int<lower=0,upper=S> site[N];  // site index
  vector[N] ppt_su;              // summer precipitation (centered) 
  vector[N] ppt_fa;              // fall precipitation   (centered)
  int<lower=0,upper=1> I0_Z;     // estimate b0_Z (1) or fix to zero (0)?
  int<lower=0,upper=1> I_su;     // estimate b_su (1) or fix to zero (0)?
  int<lower=0,upper=1> I_su_Z;   // estimate b_su_Z (1) or fix to zero (0)?
  int<lower=0,upper=1> I_fa;     // estimate b_fa (1) or fix to zero (0)?
  int<lower=0,upper=1> I_fa_Z;   // estimate b_fa_Z (1) or fix to zero (0)?
  int<lower=0> n[N];             // number of females sampled
  int<lower=0> n_psm[N];         // number of PSM females observed
  int<lower=0,upper=1> I_fit[N]; // use PSM obs i to fit model (1) or not (0)?
  int<lower=0,upper=1> I_lpd[N]; // evaluate log post predictive density for obs i (1) or not (0)?
}

transformed data {
  int<lower=1> D;                     // dimension of all site-level covariates
  matrix[S,D_normal] X_normal;        // normally distributed site-level covariates
  matrix[S,D_gamma] X_gamma;          // gamma-distributed site-level covariates
  int<lower=1,upper=N> which_fit[sum(I_fit)]; // which(I_fit == 1)
  int<lower=1,upper=N> which_lpd[sum(I_lpd)]; // which(I_lpd == 1)
  int<lower=0,upper=N> N_lpd;         // number of observations used to evaluate lpd
  
  D = D_normal + D_gamma;
  if(D_normal == 0)
    X_normal = rep_matrix(0,S,D_normal);
  else
    X_normal = X[,1:D_normal];
  if(D_gamma == 0)
    X_gamma = rep_matrix(0,S,D_gamma);
  else
    X_gamma = X[,(D_normal + 1):D];

  // Extract PSM data that will be used to fit the model
  for(i in 1:N)
  {
    if(I_fit[i] == 1)
      which_fit[sum(head(I_fit, i))] = i;
    
    if(I_lpd[i] == 1)
      which_lpd[sum(head(I_lpd, i))] = i;
  }
  
  N_lpd = sum(I_lpd);
}

parameters {
  row_vector[D] a0;           // intercepts of mu_X on link scale
  row_vector<lower=0,upper=10>[D] phi; // scale parameters of X
  vector[D*L - L*(L-1)/2] A_nid_vec;  // factor loadings (nonidentified)
  matrix[S,L] Z_nid;          // latent factor scores (nonidentified)
  real mu_b0;                 // hyper-mean of site-level PSM regression intercept
  vector[L] b0_Z_nid;         // effects of Z on site-level PSM regression intercept (nonid)
  real<lower=0,upper=10> sigma_b0; // hyper-SD of site-level PSM regression intercept
  vector[S] b0_std;           // site-level PSM regression intercepts (centered and standardized)
  real mu_b_su;               // hyper-mean of site-level summer ppt slope
  vector[L] b_su_Z_nid;       // effects of Z on site-level summer ppt slope (nonid)
  real<lower=0,upper=10> sigma_b_su; // hyper-SD of site-level summer ppt slope
  vector[S] b_su_std;         // site-level summer ppt slopes (centered and standardized)
  real mu_b_fa;               // hyper-mean of site-level fall ppt slope
  vector[L] b_fa_Z_nid;       // effects of Z on site-level fall ppt slope (nonid)
  real<lower=0,upper=10> sigma_b_fa; // hyper-SD of site-level fall ppt slope
  vector[S] b_fa_std;         // site-level fall ppt slopes (centered and standardized)
  real<lower=0,upper=10> sigma_psm; // SD of observation-level residual errors in logit(P(PSM))
  vector[N] logit_p_psm_std;  // observation-level residual errors (standardized)
}

transformed parameters {
  matrix[D,L] A_nid;             // factor loading matrix (nonidentified)
  matrix[S,D] g_mu_X;            // linear predictor of X on link scale
  vector[S] b0;                  // site-level PSM regression intercepts
  vector[S] b_su;                // site-level summer ppt slopes
  vector[S] b_fa;                // site-level fall ppt slopes
  vector[N] logit_p_psm_hat;     // expected logit P(PSM) w/o obs-level random residuals
  vector[N] logit_p_psm;         // logit P(PSM) including obs-level random residuals
  
  // Fill in factor loading matrix s.t. A_nid[1:L,1:L] is lower triangular but unconstrained
  // (note use of local variable as a counter)
  A_nid = rep_matrix(0,D,L);
  
  {
    int k;
    k = 1;
    for(i in 1:D)
      for(j in 1:min(i,L))
      {
        A_nid[i,j] = A_nid_vec[k];
        k = k + 1;
      }
  }
  
  // Calculate linear predictor matrix of X on link scale
  g_mu_X = rep_matrix(a0,S) + Z_nid*transpose(A_nid);
  
  // Shift and rescale site-level random regression coefs
  // Indicator variables (read in as data) control which coefs are
  // estimated (I = 1) or fixed to zero (I = 0)
  b0 = mu_b0 + I0_Z*Z_nid*b0_Z_nid + sigma_b0*b0_std;
  b_su = I_su*(mu_b_su + I_su_Z*Z_nid*b_su_Z_nid + sigma_b_su*b_su_std);
  b_fa = I_fa*(mu_b_fa + I_fa_Z*Z_nid*b_fa_Z_nid + sigma_b_fa*b_fa_std);
  
  // Calculate predicted P(PSM), including observation-level random residuals
  logit_p_psm_hat = b0[site] + b_su[site] .* ppt_su + b_fa[site] .* ppt_fa;
  logit_p_psm = logit_p_psm_hat + sigma_psm*logit_p_psm_std;
}

model {
  // Priors
  a0 ~ normal(0,10);
  A_nid_vec ~ normal(0,10);
  to_vector(Z_nid) ~ normal(0,1);
  mu_b0 ~ normal(0,10);
  b0_Z_nid ~ normal(0,10);
  mu_b_su ~ normal(0,10);
  b_su_Z_nid ~ normal(0,10);
  mu_b_fa ~ normal(0,10);
  b_fa_Z_nid ~ normal(0,10);
  
  
  // Observation-level residual errors in P(PSM) to allow overdispersion
  logit_p_psm_std ~ normal(0,1);
  
  // Hierarchical priors for site-level regression coefs (non-centered parameterization)
  b0_std ~ normal(0,1);   // b0 ~ normal(mu_b0 + Z_nid*b0_Z_nid, sigma_b0)
  b_su_std ~ normal(0,1); // b_su ~ normal(mu_b_su + Z_nid*b_su_Z_nid, sigma_b_su)
  b_fa_std ~ normal(0,1); // b_fa ~ normal(mu_b_fa + Z_nid*b_fa_Z_nid, sigma_b_fa)
  
  // Factor analysis likelihood of site-level covariates
  // normally distributed
  {
    matrix[S,D_normal] g_mu_X_normal;
    matrix[S,D_normal] phi_normal;
    
    if(D_normal == 0)
    {
      g_mu_X_normal = rep_matrix(0,S,D_normal);
      phi_normal = rep_matrix(0,S,D_normal);
    }
    else
    {
      g_mu_X_normal = g_mu_X[,1:D_normal];
      phi_normal = rep_matrix(head(phi, D_normal), S);
    }
    
    to_vector(X_normal) ~ normal(to_vector(g_mu_X_normal), to_vector(phi_normal));
  }
  
  // gamma-distributed
  {
    matrix[S,D_gamma] inv_mu_X_gamma;
    matrix[S,D_gamma] phi_gamma;
    
    if(D_gamma == 0)
    {
      inv_mu_X_gamma = rep_matrix(0,S,D_gamma);
      phi_gamma = rep_matrix(0,S,D_gamma);
    }
    else
    {
      inv_mu_X_gamma = exp(-g_mu_X[,(D_normal + 1):D]);
      phi_gamma = rep_matrix(tail(phi, D_gamma), S);
    }
    
    to_vector(X_gamma) ~ gamma(to_vector(phi_gamma), to_vector(phi_gamma .* inv_mu_X_gamma));
  }

  // Likelihood of observed PSM conditional on obs-level random errors
  n_psm[which_fit] ~ binomial_logit(n[which_fit], logit_p_psm[which_fit]);
}

generated quantities {
  vector[L] sign_A_diag;  // sign of upper 1:L diagonal of nonidentified loading matrix
  matrix[D,L] A;          // loading matrix (identified)
  matrix[S,L] Z;          // factor scores (identified)
  vector[L] b0_Z;         // effects of Z on site-level PSM regression intercept (identified)
  vector[L] b_su_Z;       // effects of Z on site-level summer precip slope (identified)
  vector[L] b_fa_Z;       // effects of Z on site-level fall precip slope (identified)
  vector[N] p_psm;        // predicted P(PSM)
  int N_MC;               // sample size for Monte Carlo estimate of marginal likelihood
  vector[N_lpd] ll_psm;   // pointwise marginal log-likelihood of PSM obs

  // Reflect factors and loadings to ensure identifiability
  // Identified loading matrix A has A[1:L,1:L] lower triangular with positive diagonal
  for(j in 1:L)
    sign_A_diag[j] = A_nid[j,j] > 0 ? 1 : -1;
  
  A = diag_post_multiply(A_nid, sign_A_diag);
  Z = diag_post_multiply(Z_nid, sign_A_diag);
  b0_Z = b0_Z_nid .* sign_A_diag;
  b_su_Z = b_su_Z_nid .* sign_A_diag;
  b_fa_Z = b_fa_Z_nid .* sign_A_diag;
  
  // Transform P(PSM) from logit scale
  p_psm = 1 ./ (1 + exp(-logit_p_psm));

  // Pointwise marginal log-likelihood of PSM observations 
  // for use in model comparison via LOO or WAIC
  // Marginal likelihood is approximated by Monte Carlo integration over the
  // observation-level random residuals
  N_MC = 1000;  // hard-coded
  ll_psm = rep_vector(0,N_lpd);
  for(i in 1:N_lpd)  // will only execute if N_lpd > 0
  {
    vector[N_MC] ll_psm_MC;
    
    for(j in 1:N_MC)
    {
      real logit_p_psm_MC;
      
      logit_p_psm_MC = normal_rng(logit_p_psm_hat[which_lpd[i]], sigma_psm);
      ll_psm_MC[j] = binomial_logit_lpmf(n_psm[which_lpd[i]] | n[which_lpd[i]], logit_p_psm_MC);
    }
    
    ll_psm[i] = log_sum_exp(ll_psm_MC) - log(N_MC);
  }
}
