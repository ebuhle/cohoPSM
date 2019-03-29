model
{
  #-------------------------------------------------------------------
  # PRIORS
  #-------------------------------------------------------------------
  
  for(i in 1:D)
  {
    a0[i] ~ dnorm(0,1e-6)   # intercept of X[,i]
    
    # DxL matrix of factor loadings constrained s.t. A[1:L,1:L] is lower triangular 
    # with positive diagonal.
    for(j in 1:L)
    {
      AA[i,j] ~ dnorm(0,1e-6)    
      A.tilde[i,j] <- AA[i,j]*ifelse(j > i, 0, 1) 
      A[i,j] <- A.tilde[i,j]*ifelse(A.tilde[j,j] < 0, -1, 1)*sigma.Z[j]  # transform to identified loadings
      #       A[i,j] <- AA[i,j]*ifelse(j > i, 0, 1)*ifelse(AA[j,j] < 0, -1, 1)
    }
    
    phi[i] ~ dunif(0,100)  # shape (a) parameter of X[,i]
  }
  
  # SD of Z.tilde (redundant parameter-expanded version of factors) 
  for(j in 1:L)
  {
    tau.Z[j] ~ dgamma(0.01,0.01)
    sigma.Z[j] <- 1/sqrt(tau.Z[j])
  }
  
  #SxL matrix of latent factor scores.
  for(i in 1:S)
  {
    for(j in 1:L)
    {
      #       Z.tilde[i,j] ~ dnorm(0,1)
      #       Z[i,j] <- Z.tilde[i,j]*ifelse(AA[j,j] < 0, -1, 1)
      Z.tilde[i,j] ~ dnorm(0,tau.Z[j])        
      Z[i,j] <- Z.tilde[i,j]*ifelse(A.tilde[j,j] < 0, -1, 1)/sigma.Z[j]  # transform to identfied factor scores
    }
  }
  
  #-------------------------------------------------------------------
  # LATENT-FACTOR MODEL FOR COVARIATES
  #-------------------------------------------------------------------
  
  for(i in 1:S)
  {
    for(j in 1:D)
    {
      mu[i,j] <- exp(a0[j] + sum(A.tilde[j,1:L]*Z.tilde[i,1:L]))
      #       mu[i,j] <- exp(a0[j] + sum(A.tilde[j,1:L]*Z[i,1:L]))
      
      # Likelihood of observed X[i,j]
      X[i,j] ~ dgamma(phi[j], phi[j]/mu[i,j])
    }
  }
}
