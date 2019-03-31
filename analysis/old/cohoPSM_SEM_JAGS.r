#####################################################################################################
# MULTILEVEL LOGISTIC MODEL FOR COHO PSM USING SUPERVISED DIMENSION REDUCTION VIA LATENT FACTORS
#####################################################################################################


model
{
  #-------------------------------------------------------------------
  # PRIORS
  #-------------------------------------------------------------------
  
  for(i in 1:D)
  {
    a0[i] ~ dnorm(0,1e-6)  # intercept of X[,i] on link scale
    phi[i] ~ dunif(0,1e6)  # shape parameter of X[,i]
  }
  
  for(j in 1:L)
  {
    # DxL matrix of factor loadings constrained s.t. A[1:L,1:L] is 
    # lower triangular with positive diagonal.
    for(i in 1:D)
    {
      AA[i,j] ~ dnorm(0,1e-6)                                      # unconstrained loadings
      A.px[i,j] <- AA[i,j]*ifelse(j > i, 0, 1)                     # parameter-expanded loadings
      A[i,j] <- A.px[i,j]*ifelse(A.px[j,j] < 0, -1, 1)*sigma.Z[j]  # identified loadings
    }
    
    # precision of Z.px (redundant parameter-expanded version of latent factors) 
    tau.Z[j] ~ dgamma(0.01,0.01)
    sigma.Z[j] <- 1/sqrt(tau.Z[j])

    # SxL matrix of latent factor scores
    for(i in 1:S)
    {
      Z.px[i,j] ~ dnorm(0,tau.Z[j])                                # parameter-expanded factor scores      
      Z[i,j] <- Z.px[i,j]*ifelse(A.px[j,j] < 0, -1, 1)/sigma.Z[j]  # identfied factor scores
    }
    
    B0.Z.px[j] ~ dnorm(0,1e-6)                                     # effect of Z on site-level PSM intercept (PX version)
    B0.Z[j] <- B0.Z.px[j]*ifelse(A.px[j,j] < 0, -1, 1)*sigma.Z[j]  # identified version
    
    B.su.Z.px[j] ~ dnorm(0,1e-6)                                   # effect of Z on site-level su ppt slope (PX)
    B.su.Z[j] <- B.su.Z.px[j]*ifelse(A.px[j,j] < 0, -1, 1)*sigma.Z[j] # identified version
    
    B.fa.Z.px[j] ~ dnorm(0,1e-6)                                   # effect of Z on site-level fa ppt slope (PX)
    B.fa.Z[j] <- B.fa.Z.px[j]*ifelse(A.px[j,j] < 0, -1, 1)*sigma.Z[j] # identified version
  }
    
  B0 ~ dnorm(0,1e-6)         # hyper-mean of PSM regression intercept
  B.su ~ dnorm(0,1e-6)       # hyper-mean of summer ppt slope (i.e., when Z = 0)
  B.fa ~ dnorm(0,1e-6)       # hyper-mean of fall ppt slope (i.e., when Z = 0)
  sigma.b0 ~ dunif(0,10)     # hyper-SD of site-level PSM intercepts
  sigma.b.su ~ dunif(0,10)   # hyper-SD of site-level summer ppt slopes
  sigma.b.fa ~ dunif(0,10)   # hyper-SD of site-level fall ppt slopes  
  
  # Site-level coefficients for multilevel PSM regression
  for(i in 1:S)
  {
    b0[i] ~ dnorm(B0 + sum(B0.Z.px*Z.px[i,1:L]), pow(sigma.b0,-2))          # site i PSM intercept
    b.su[i] ~ dnorm(B.su + sum(B.su.Z.px*Z.px[i,1:L]), pow(sigma.b.su,-2))  # site i summer ppt slope
    b.fa[i] ~ dnorm(B.fa + sum(B.fa.Z.px*Z.px[i,1:L]), pow(sigma.b.fa,-2))  # site i fall ppt slope
  }
    
  #-------------------------------------------------------------------
  # LATENT-FACTOR MODEL FOR COVARIATES
  #-------------------------------------------------------------------

  for(i in 1:S)
  {
    # Normally distributed covariates
    for(j in normal.indx)
    {
      mu[i,j] <- a0[j] + sum(A.px[j,1:L]*Z.px[i,1:L])
      X[i,j] ~ dnorm(mu[i,j], pow(phi[j],-2))
    }

    # Gamma-distributed covariates
    for(j in gamma.indx)
    {
      mu[i,j] <- exp(a0[j] + sum(A.px[j,1:L]*Z.px[i,1:L]))
      X[i,j] ~ dgamma(phi[j], phi[j]/mu[i,j])
    }
  }
  
  #-------------------------------------------------------------------
  # MULTILEVEL LOGISTIC REGRESSION MODEL FOR PSM
  #-------------------------------------------------------------------
  
  for(i in 1:length(n))
  {
    p.psm[i] <- ilogit(b0[site[i]] + b.su[site[i]]*ppt.su[i] + b.fa[site[i]]*ppt.fa[i])
    n.psm[i] ~ dbinom(min(max(p.psm[i], 1e-6), 1 - 1e-6), n[i])
  }
}
