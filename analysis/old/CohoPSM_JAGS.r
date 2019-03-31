#======================================================
# COHO PRESPAWN MORTALITY ANALYSIS
#
# Mixed-effects regression models using hierearchical
# Bayesian analysis, with Bayesian model selection
# (Bayes factors).
#======================================================


#========================================================================================================================================
# ANALYSIS WITH RESTRICTED LAND USE / LAND COVER PREDICTORS AND PRECIPITATION
#========================================================================================================================================

# Response is binomial (n.psm/n.total)
# There is a random effect on the intercept, grouped by site
# Candidate models consist of combinations of the predictors:
# Impervious, FC00.Road, dry, rain, Impervious:dry, Impervious:rain, FC00.Road:dry, FC00.Road:rain
# with the restriction that Impervious OR FC00.Road may appear in a model, but not both.
# Premise: compare hypotheses that each of these factors is the best predictor of PSM.

library(R2jags)
library(lme4)

# Data objects
N <- nrow(psm)
S <- length(levels(psm$site))
site <- as.numeric(psm$site)
Impervious <- (psm$Impervious - mean(psm$Impervious))/sd(psm$Impervious)
FC00.Road <- (psm$FC00.Road - mean(psm$FC00.Road))/sd(psm$FC00.Road)
dry <- (psm$dry - mean(psm$dry))/sd(psm$dry)
rain <- (psm$rain - mean(psm$rain))/sd(psm$rain)
n.total <- psm$n.total
n.psm <- psm$n.psm


# Write out model file
mod.psm <-
"model
{
	# Site-specific intercepts
	for(j in 1:S)
	{
	b0[j] ~ dnorm(mu.b0, tau.b0)
	}

	# Fitted values and binomial likelihood
	for(i in 1:N)
	{
	fit.logit.psm[i] <- b0[site[i]] + w.Imp*b.Imp*Impervious[i] + w.Road*b.Road*FC00.Road[i] + w.dry*b.dry*dry[i] + w.rain*b.rain*rain[i] + 
				w.Impdry*b.Impdry*Impervious[i]*dry[i] + w.Imprain*Impervious[i]*rain[i] +
				w.Roaddry*FC00.Road[i]*dry[i] + w.Roadrain*FC00.Road[i]*rain[i]
	fit.psm[i] <- 1/(1 + exp(-fit.logit.psm[i]))
	n.psm[i] ~ dbin(fit.psm[i], n.total[i])
	}

	# Priors
	mu.b0 ~ dnorm(0,1.0E-6)
	sigma.b0 ~ dunif(0,100)
	tau.b0 <- pow(sigma.b0,-2)
	w.Imp ~ dbern(0.5)
	b.Imp ~ dnorm(0,1.0E-6)
	w.Road ~ dbern(0.5)
	b.Road ~ dnorm(0,1.0E-6)
	w.dry ~ dbern(0.5)
	b.dry ~ dnorm(0,1.0E-6)
	w.rain ~ dbern(0.5)
	b.rain ~ dnorm(0,1.0E-6)
	v.Impdry <- 0 #~ dbern(0.5)
	w.Impdry <- w.Imp*w.dry*v.Impdry
	b.Impdry <- 0 #~ dnorm(0,1.0E-6)
	v.Imprain <- 0 #~ dbern(0.5)
	w.Imprain <- w.Imp*w.rain*v.Imprain
	b.Imprain <- 0 #~ dnorm(0,1.0E-6)
	v.Roaddry <- 0 #~ dbern(0.5)
	w.Roaddry <- w.Road*w.dry*v.Roaddry
	b.Roaddry <- 0 #~ dnorm(0,1.0E-6)
	v.Roadrain <- 0 #~ dbern(0.5)
	w.Roadrain <- w.Road*w.rain*v.Roadrain
	b.Roadrain <- 0 #~ dnorm(0,1.0E-6)
}"
write(mod.psm, "mod.psm.txt")
rm(mod.psm)

# Fit model
inits.psm <- function() list(mu.b0 = rnorm(1,0,1), sigma.b0 = runif(1,0,1), w.Imp = 1, b.Imp = rnorm(1,0,1),
w.Road = 1, b.Road = rnorm(1,0,1), w.dry = 1, b.dry = rnorm(1,0,1), 
w.rain = 1, b.rain = rnorm(1,0,1), 
v.Impdry = 1, b.Impdry = rnorm(1,0,1), v.Imprain = 1, b.Imprain = rnorm(1,0,1), 
v.Roaddry = 1, b.Roaddry = rnorm(1,0,1), v.Roadrain = 1, b.Roadrain = rnorm(1,0,1))

inits.psm <- function() list(mu.b0 = rnorm(1,0,1), sigma.b0 = runif(1,0,1), w.Imp = 1, b.Imp = rnorm(1,0,1),
w.Road = 1, b.Road = rnorm(1,0,1), w.dry = 1, b.dry = rnorm(1,0,1), 
w.rain = 1, b.rain = rnorm(1,0,1))

jags.psm <- jags2(data = c("N","S","site","Impervious","FC00.Road","dry","rain","n.total","n.psm"), inits = inits.psm, 
	parameters.to.save = c(names(inits.psm()), "b0"), #, "w.Impdry", "w.Imprain", "w.Roaddry", "w.Roadrain"), 
	model.file = "mod.psm.txt", n.chains = 2, n.iter = 20000, n.burnin = 2000, n.thin = 100, DIC = T, refresh = 100)

# Basic plots
print(jags.psm,2);plot(jags.psm)
plot(mcmc.list(mcmc(jags.psm$sims.array[,1,]),mcmc(jags.psm$sims.array[,2,])),ask=T)
autocorr.plot(mcmc.list(mcmc(jags.psm$sims.array[,1,]),mcmc(jags.psm$sims.array[,2,])),ask=T)








