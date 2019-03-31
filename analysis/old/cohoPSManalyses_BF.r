library(jagsUI)
library(coda)


#==================================================================
# DATA
#==================================================================

#----------------------------------------------
# PSM SAMPLING SITES FOR ANALYSIS
#----------------------------------------------

# read in spawner data
spawner.data <- read.csv("spawner.data.csv", header=T)

# read in landscape data
spatial.data <- read.csv("spatial.data.csv", header=T)

# LU/LC and roads data
lulc.roads.data <- spatial.data[,c("site","area","ccap.decid","ccap.open","ccap.evgrn","ccap.hdev","ccap.ldev",
                                   "ccap.mdev","ccap.mxforest","ccap.wetland","ccap.ag","ccap.other",
                                   "roads.1","roads.2","roads.3","roads.4","roads.5","traffic","nlcd.imperv",
                                   "restoration","pop.census","pop.lscan")]

# convert basin area from m2 to ha
lulc.roads.data$area <- lulc.roads.data$area/10e4

# convert to km of road/km2
lulc.roads.data[,substring(names(lulc.roads.data),1,5)=="roads"] <- lulc.roads.data[,substring(names(lulc.roads.data),1,5)=="roads"]/1000

# matrix of annual summer and fall precip with sites and years as rows
ppt.data <- spatial.data[,substring(names(spatial.data),1,8) %in% c("site","ppt.su.2","ppt.fa.2")]
ppt.su.cols <- substring(names(ppt.data),1,8)=="ppt.su.2"
ppt.fa.cols <- substring(names(ppt.data),1,8)=="ppt.fa.2"
n.yrs.ppt <- sum(ppt.su.cols)

ppt.data <- data.frame(site=rep(ppt.data$site, each=n.yrs.ppt), 
                       year=rep(substring(names(ppt.data)[ppt.su.cols], 8), nrow(ppt.data)),
                       ppt.su=as.vector(t(as.matrix(ppt.data[,ppt.su.cols]))),
                       ppt.fa=as.vector(t(as.matrix(ppt.data[,ppt.fa.cols]))))

# merge basin-specific variables and annual su and fa ppt into spawner data
psm <- merge(spawner.data, ppt.data, by=c("site","year"))
spatial.data.cols <- substring(names(spatial.data),1,3) != "ppt"
psm <- merge(psm, spatial.data[,spatial.data.cols], by=c("ID","site"))
psm <- data.frame(psm[,c("ID","site","watershed","year","n","n.psm","ppt.su","ppt.fa","area")], psm[,-(1:10)])
psm <- psm[order(psm$site,psm$year),]

#----------------------------------------------
# ADDITIONAL SITES FOR PSM PREDICTIONS
#----------------------------------------------

# read in landscape data
spatial.data.pre <- read.csv("spatial.data.predict.csv", header=T)
names(spatial.data.pre)[names(spatial.data.pre)=="ID"] <- "site"

# convert basin area from m2 to ha
spatial.data.pre$area <- spatial.data.pre$area/10e4

# convert to km of road/km2
spatial.data.pre[,substring(names(spatial.data.pre),1,5)=="roads"] <- spatial.data.pre[,substring(names(spatial.data.pre),1,5)=="roads"]/1000

spatial.data.pre <- spatial.data.pre[spatial.data.pre$coho==1,] # only include basins w/ coho 
spatial.data.pre$site <- factor(spatial.data.pre$site)

lulc.roads.data.pre <- spatial.data.pre[,names(lulc.roads.data)]

# Combine PSM-sampling and unsampled-site data frames
lulc.roads.data.all <- rbind(data.frame(data="psm", lulc.roads.data), 
                             data.frame(data="pre", lulc.roads.data.pre))

# Add rows to psm for unobserved sites
psm.all <- as.data.frame(matrix(nrow=nrow(spatial.data.pre), ncol=ncol(psm)))
names(psm.all) <- names(psm)
psm.all$site <- lulc.roads.data.pre$site
psm.all$n <- 1
psm.all$ppt.su <- mean(psm$ppt.su)
psm.all$ppt.fa <- mean(psm$ppt.fa)
psm.all <- rbind(data.frame(data="psm", psm), 
                 data.frame(data="pre", psm.all))


#==================================================================
# LATENT-VARIABLE REGRESSION MODELS FOR PSM
#==================================================================

#------------------------------------------------------
# Fit model to sites with PSM observations
#------------------------------------------------------

# # site-level covariates
# X <- as.matrix(lulc.roads.data[,c("ccap.hdev","ccap.ldev","ccap.mdev",
#                                   "ccap.decid","ccap.evgrn","ccap.mxforest",
#                                   "ccap.open","ccap.wetland","ccap.ag", "ccap.other",
#                                   "roads.1","roads.2","roads.3","roads.4","roads.5","traffic","nlcd.imperv",
#                                   "restoration","pop.census","pop.lscan")])
# 
# # indices of covariates assumed to be beta (gamma) distributed
# beta.indx <- which(substring(dimnames(X)[[2]],1,4) == "nlcd")
# dirch.indx <- which(substring(dimnames(X)[[2]],1,4) == "ccap")
# # beta.indx <- which(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd"))
# gamma.indx <- which(!(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd")))
# 
# # rescale gamma-distributed vars to SD = 1, 
# # bound all vars away from 0 and beta- and Dirichlet-distributed vars away from 1,
# # and force Dirichlet vars to sum to 1
# X[,gamma.indx] <- sweep(X[,gamma.indx], 2, apply(X[,gamma.indx], 2, sd), "/")
# X[X==0] <- 1e-4
# # X[,beta.indx][X[,beta.indx]==1] <- 1 - 1e-4
# X[,c(beta.indx,dirch.indx)][X[,c(beta.indx,dirch.indx)]==1] <- 1 - 1e-4
# X[,dirch.indx] <- sweep(X[,dirch.indx], 1, rowSums(X[,dirch.indx]), "/")
# 
# S <- nrow(X)                     # number of sites
# D <- ncol(X)                     # dimension of site-level covariates
# L <- 1                           # dimension of latent-factor space (user-specified!)
# n <- psm$n                       # number of females sampled
# n.psm <- psm$n.psm               # number of PSM females
# site <- as.numeric(psm$site)     # site index
# 
# ppt.su <- as.vector(scale(psm$ppt.su/10, scale=F))    # center and convert mm to cm
# ppt.fa <- as.vector(scale(psm$ppt.fa/10, scale=F))

# function to generate initial values for chains
inits.fn <- function() list(a0 = rnorm(D, sapply(1:D, function(j) ifelse(j %in% beta.indx,
                                                                         qlogis(mean(X[,j])),
                                                                         log(mean(X[,j])))), 1),
                            AA = matrix(rnorm(D*L, 0, 1), D, L),
                            phi = sapply(1:D, function(j) ifelse(j %in% beta.indx, 
                                                                 runif(1, 5, 10), 
                                                                 runif(1, 0.5, 2))),
                            B0 = rnorm(1, qlogis(mean(n.psm/n, na.rm=T)), 1),
                            B0.Z.px = rnorm(L, 0, 0.1),
                            sigma.b0 = runif(1, 0.1, 1),
                            B.su = rnorm(1, 0, 1),
                            B.su.Z.px = rnorm(L, 0, 0.1),
                            sigma.b.su = runif(1, 0.1, 1),
                            B.fa = rnorm(1, 0, 1),
                            B.fa.Z.px = rnorm(L, 0, 0.1),
                            sigma.b.fa = runif(1, 0.1, 1))

# # Fit it!
# jags.psm <- jags(data = c("beta.indx","dirch.indx","gamma.indx","X","S","D","L","n","n.psm","site","ppt.su","ppt.fa"), 
#                  inits = inits.fn, 
#                  parameters.to.save = c("a0","A","Z","phi","mu",
#                                         "B0","B0.Z","sigma.b0","b0",
#                                         "B.su","B.su.Z","sigma.b.su","b.su",
#                                         "B.fa","B.fa.Z","sigma.b.fa","b.fa",
#                                         "p.psm"), 
#                  model.file = "cohoPSM_SEM_JAGS.r",
#                  n.chains = 3, n.iter = 100000, n.burnin = 20000, n.thin = 80, 
#                  DIC = T, parallel = TRUE, seed = Sys.time())
# 
# 
# # Plot and print summaries of results
# print(jags.psm,2)
# 
# # Plot chains and histograms
# plot(jags.psm)


#---------------------------------------------------------
# Fit model to sites with PSM observations plus
# unsampled sites to generate predictions for the latter
#---------------------------------------------------------

# site-level covariates
X <- as.matrix(lulc.roads.data.all[,c("ccap.hdev","ccap.ldev","ccap.mdev",
                                      "ccap.decid","ccap.evgrn","ccap.mxforest",
                                      "ccap.open","ccap.wetland","ccap.ag", "ccap.other",
                                      "roads.1","roads.2","roads.3","roads.4","roads.5","traffic","nlcd.imperv",
                                      "restoration","pop.census","pop.lscan")])

# indices of covariates assumed to be beta (gamma) distributed
beta.indx <- which(substring(dimnames(X)[[2]],1,4) == "nlcd")
dirch.indx <- which(substring(dimnames(X)[[2]],1,4) == "ccap")
# beta.indx <- which(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd"))
gamma.indx <- which(!(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd")))

# rescale gamma-distributed vars to SD = 1, 
# bound all vars away from 0 and beta- and Dirichlet-distributed vars away from 1,
# and force Dirichlet vars to sum to 1
X[,gamma.indx] <- sweep(X[,gamma.indx], 2, apply(X[,gamma.indx], 2, sd), "/")
X[X==0] <- 1e-4
# X[,beta.indx][X[,beta.indx]==1] <- 1 - 1e-4
X[,c(beta.indx,dirch.indx)][X[,c(beta.indx,dirch.indx)]==1] <- 1 - 1e-4
X[,dirch.indx] <- sweep(X[,dirch.indx], 1, rowSums(X[,dirch.indx]), "/")

S <- nrow(X)                     # number of sites
D <- ncol(X)                     # dimension of site-level covariates
L <- 1                           # dimension of latent-factor space (user-specified!)
n <- psm.all$n                       # number of females sampled
n.psm <- psm.all$n.psm               # number of PSM females
site <- as.numeric(psm.all$site)     # site index

ppt.su <- as.vector(scale(psm.all$ppt.su/10, scale=F))    # center and convert mm to cm
ppt.fa <- as.vector(scale(psm.all$ppt.fa/10, scale=F))

# Fit it, using jags.basic() to return a minimal mcmc.list object
jags.psm.all <- jags.basic(data = c("beta.indx","dirch.indx","gamma.indx","X","S","D","L","n","n.psm","site","ppt.su","ppt.fa"), 
                           inits = inits.fn, 
                           parameters.to.save = c("a0","A","Z","phi","mu",
                                                  "B0","B0.Z","sigma.b0","b0",
                                                  "B.su","B.su.Z","sigma.b.su","b.su",
                                                  "B.fa","B.fa.Z","sigma.b.fa","b.fa",
                                                  "p.psm"), 
                           model.file = "cohoPSM_SEM_JAGS.r",
                           #                            n.chains = 3, n.iter = 100000, n.burnin = 20000, n.thin = 80, 
                           n.chains = 3, n.iter = 200, n.burnin = 100, n.thin = 1, 
                           DIC = T, parallel = TRUE, seed = Sys.time(), save.model = FALSE)

# # Plot chains and histograms (if you dare)
# plot(jags.psm.all)
# plot(jags.psm.all[,substring(varnames(jags.psm.all),1,4)=="b.su"], ask=T)

# Store predicted P(PSM) in matrix
psm.pre <- as.matrix(jags.psm.all[,substring(varnames(jags.psm.all),1,5)=="p.psm"])
row.names(psm.pre) <- NULL
psm.pre <- data.frame(site=psm.all$site[psm.all$data=="pre"], 
                      psm.mean=colMeans(psm.pre)[psm.all$data=="pre"],
                      logit.psm.sd=apply(qlogis(psm.pre),2,sd)[psm.all$data=="pre"])
write.table(psm.pre, "PSM_predictions.txt", sep="\t", row.names=FALSE)

# Save workspace
save.image(file="cohoPSMpredict.RData")