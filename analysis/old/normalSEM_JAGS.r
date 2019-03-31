model
{
  #-------------------------------------------------------------------
  # PRIORS
  #-------------------------------------------------------------------
  
  for(i in 1:D)
  {
    a0[i] ~ dnorm(0,1e-6)         # intercept of X[,i]
    #     tau.X[i] ~ dgamma(0.001,0.001)  # precision of X[,i]
    #     phi[i] <- 1/sqrt(tau.X[i])      # SD
    sigma.X[i] ~ dunif(0,10)
    tau.X[i] <- pow(sigma.X[i],-2)
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
    
    # SxL matrix of latent factor scores.
    for(i in 1:S)
    {
      Z.px[i,j] ~ dnorm(0,tau.Z[j])                                # parameter-expanded factor scores      
      Z[i,j] <- Z.px[i,j]*ifelse(A.px[j,j] < 0, -1, 1)/sigma.Z[j]  # identfied factor scores
    }
  }

  #-----------------------------------------------------------------------
  # LATENT-FACTOR MODEL FOR COVARIATES AND REGRESSION MODEL FOR RESPONSE
  #-----------------------------------------------------------------------
  
  for(i in 1:length(y))
  {
    for(j in 1:D)
    {
      mu.X[i,j] <- a0[j] + sum(A.px[j,1:L]*Z.px[i,1:L])  # mean of covariate j
      X[i,j] ~ dnorm(mu[i,j], tau.X[j])                  # likelihood of covariate j
    }
    
    mu.y[i] <- b0 + sum(b.px[1:L]*Z.px[i,1:L])           # mean of y
    y[i] ~ dnorm(mu.y[i], tau.y)
  }
}
