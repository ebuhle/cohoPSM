model
{
  #-------------------------------------------------------------------
  # PRIORS
  #-------------------------------------------------------------------
  
  for(i in 1:D)
  {
    a0[i] ~ dnorm(0,1e-6)   # intercept of X[,i] on link scale
    
    # DxL matrix of factor loadings constrained s.t. A[1:L,1:L] is lower triangular 
    # with positive diagonal.
    for(j in 1:L)
    {
      AA[i,j] ~ dnorm(0,1e-6)    
      A.tilde[i,j] <- AA[i,j]*ifelse(j > i, 0, 1) 
      A[i,j] <- A.tilde[i,j]*ifelse(A.tilde[j,j] < 0, -1, 1)*sigma.Z[j]  # transform to identified loadings
    }
    
    phi[i] ~ dlnorm(0,1/25)  # dispersion parameter of X[,i]
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
      Z.tilde[i,j] ~ dnorm(0,tau.Z[j])        
      Z[i,j] <- Z.tilde[i,j]*ifelse(A.tilde[j,j] < 0, -1, 1)/sigma.Z[j]  # transform to identfied factor scores
    }
  }
  
  #-------------------------------------------------------------------
  # LATENT-FACTOR MODEL FOR COVARIATES
  #-------------------------------------------------------------------
  
  for(i in 1:S)
  {
    for(j in 1:length(beta.indx))
    {
      mu[i,beta.indx[j]] <- ilogit(a0[beta.indx[j]] + sum(A.tilde[beta.indx[j],1:L]*Z.tilde[i,1:L]))
      
      # Likelihood of observed X[i,j]
      X[i,beta.indx[j]] ~ dbeta(mu[i,beta.indx[j]]*phi[beta.indx[j]], (1 - mu[i,beta.indx[j]])*phi[beta.indx[j]])
    }
    
    for(j in 1:length(gamma.indx))
    {
      mu[i,gamma.indx[j]] <- exp(a0[gamma.indx[j]] + sum(A.tilde[gamma.indx[j],1:L]*Z.tilde[i,1:L]))
      
      # Likelihood of observed X[i,j]
      X[i,gamma.indx[j]] ~ dgamma(phi[gamma.indx[j]], phi[gamma.indx[j]]/mu[i,gamma.indx[j]])
    }
  }
}
