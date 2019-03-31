S <- 50
D <- 10
L <- 1

Z.true <- matrix(rnorm(S*L), S, L)
a0.true <- rnorm(D)
A.true <- matrix(c(abs(runif(1,-2,2)), runif(D*L - 1, -2, 2)), D, L)
# mu.true <- matrix(a0.true, S, D, byrow=T) + Z.true %*% t(A.true)
phi.true <- runif(D, 0.5, 2)
mu.true <- matrix(NA,S,D)
X <- matrix(NA,S,D)
for(i in 1:S)
{
  for(j in 1:D)
  {
    mu.true[i,j] <- exp(a0.true[j] + sum(A.true[j,1:L]*Z.true[i,1:L]))
    X[i,j] <- pmax(rgamma(1, phi.true[j], phi.true[j]/mu.true[i,j]), 1e-6)
  }
}

# function to generate initial values for chains
inits.fn <- function() list(a0 = rnorm(D, colMeans(log(X)), 1),
                            AA = matrix(rnorm(D*L, 0, 1), D, L),
                            phi = runif(D, 0.5, 2))

# Fit it!
jags.FA <- jags.parallel(data = c("X","S","D","L"), 
                         inits = inits.fn, 
                         parameters.to.save = c("a0","A","Z","phi","mu"), 
                         model.file = "gammaFA_JAGS.r",
                         n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 5, 
                         DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.FA,2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x$BUGS$sims.array[,d,substring(dimnames(x$BUGS$sims.array)[[3]],1,2) != "mu"]), 
                      1:jags.FA$BUGS$n.chains, list(x=jags.FA), SIMPLIFY=F)), ask=T)

# True vs. fitted a0, A, Z, phi
par(mfrow=c(2,3))
plot(a0.true, jags.FA$BUGS$mean$a0);abline(0,1)
plot(A.true, jags.FA$BUGS$mean$A);abline(0,1)
plot(Z.true, jags.FA$BUGS$mean$Z);abline(0,1)
plot(phi.true, jags.FA$BUGS$mean$phi);abline(0,1)
plot(mu.true, jags.FA$BUGS$mean$mu);abline(0,1)


