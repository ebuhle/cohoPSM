###############################################################################################################
## HIERARCHICAL LOGISTIC REGRESSION MODEL FOR COHO PRESPAWN MORTALITY



model
{
  # PRIORS
  
  # Site-level regression coefs
  b0.0 ~ dnorm(0,1e-6)              # intercept of site-level intercept
  b0.lulc ~ dnorm(0,1e-6)           # LU/LC slope of site-level intercept
  sigma.b0 ~ dunif(0,10)            # SD of site-level intercept
  
  b.ppt.su.0 ~ dnorm(0,1e-6)        # intercept of site-level summer precip slope
  b.ppt.su.lulc ~ dnorm(0,1e-6)     # LU/LC slope of site-level summer precip slope
  sigma.b.ppt.su ~ dunif(0,10)      # SD of site-level summer precip slope
  
  b.ppt.fa.0 ~ dnorm(0,1e-6)        # intercept of site-level fall precip slope
  b.ppt.fa.lulc ~ dnorm(0,1e-6)     # LU/LC slope of site-level fall precip slope
  sigma.b.ppt.fa ~ dunif(0,10)      # SD of site-level fall precip slope
  
  b.ppt.su.fa.0 ~ dnorm(0,1e-6)     # intercept of site-level fall precip x summer precip interaction slope
  b.ppt.su.fa.lulc ~ dnorm(0,1e-6)  # LU/LC slope of site-level fall precip x summer precip interaction slope
  sigma.b.ppt.su.fa ~ dunif(0,10)   # SD of site-level fall precip x summer precip interaction slope

  # SITE-LEVEL MODEL
  
  for(j in 1:max(site))
  {
    b0[j] ~ dnorm(b0.0 + b0.lulc*lulc[j], pow(sigma.b0,-2))  # site-level intercept    
    #     b.ppt.su[j] <- b.ppt.su.0 + b.ppt.su.lulc*lulc[j]          # site-level summer precip slope
    #     b.ppt.fa[j] <- b.ppt.fa.0 + b.ppt.fa.lulc*lulc[j]          # site-level fall precip slope
    b.ppt.su[j] ~ dnorm(b.ppt.su.0 + b.ppt.su.lulc*lulc[j], pow(sigma.b.ppt.su,-2))  # site-level summer precip slope
    b.ppt.fa[j] ~ dnorm(b.ppt.fa.0 + b.ppt.fa.lulc*lulc[j], pow(sigma.b.ppt.fa,-2))  # site-level fall precip slope
    #     b.ppt.su.fa[j] ~ dnorm(b.ppt.su.fa.0, pow(sigma.b.ppt.su.fa,-2)) # site-level su x fa precip intxn slope
    b.ppt.su.fa[j] ~ dnorm(b.ppt.su.fa.0 + b.ppt.su.fa.lulc*lulc[j], pow(sigma.b.ppt.su.fa,-2)) # site-level su x fa precip intxn slope
  }
  
  # OBSERVATION-LEVEL MODEL
  
  for(i in 1:length(n))
  {
    lp.psm[i] <- b0[site[i]] + b.ppt.su[site[i]]*ppt.su[i] + b.ppt.fa[site[i]]*ppt.fa[i] + b.ppt.su.fa[site[i]]*ppt.su[i]*ppt.fa[i]  # predicted logit(P(PSM))
    p.psm[i] <- ilogit(lp.psm[i])                  # inverse logit transform
    n.psm[i] ~ dbinom(p.psm[i], n[i])              # likelihood of observed PSM count
  }
  
  # PREDICTIONS FOR UNOBSERVED SITES
  for(i in 1:length(lulc.pre))
  {
    b0.pre[i] ~ dnorm(b0.0 + b0.lulc*lulc.pre[i], pow(sigma.b0,-2))  # site-level intercept
    #     b.ppt.su.pre[i] <- b.ppt.su.0 + b.ppt.su.lulc*lulc.pre[i]          # site-level summer precip slope
    #     b.ppt.fa.pre[i] <- b.ppt.fa.0 + b.ppt.fa.lulc*lulc.pre[i]          # site-level fall precip slope
    b.ppt.su.pre[i] ~ dnorm(b.ppt.su.0 + b.ppt.su.lulc*lulc.pre[i], pow(sigma.b.ppt.su,-2))  # site-level summer precip slope
    b.ppt.fa.pre[i] ~ dnorm(b.ppt.fa.0 + b.ppt.fa.lulc*lulc.pre[i], pow(sigma.b.ppt.fa,-2))  # site-level fall precip slope
    #     b.ppt.su.fa.pre[i] ~ dnorm(b.ppt.su.fa.0, pow(sigma.b.ppt.su.fa,-2)) # site-level su x fa precip intxn slope
    b.ppt.su.fa.pre[i] ~ dnorm(b.ppt.su.fa.0 + b.ppt.su.fa.lulc*lulc.pre[i], pow(sigma.b.ppt.su.fa,-2)) # site-level su x fa precip intxn slope
    
    lp.psm.pre[i] <- b0.pre[i] + b.ppt.su.pre[i]*ppt.su.pre[i] + b.ppt.fa.pre[i]*ppt.fa.pre[i] + b.ppt.su.fa.pre[i]*ppt.su.pre[i]*ppt.fa.pre[i]  # predicted logit(P(PSM))
    p.psm.pre[i] <- ilogit(lp.psm.pre[i])          # inverse logit transform
  }
}

