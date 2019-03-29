options(device=windows)
# library(lme4)
# library(R2jags)
library(jagsUI)
library(coda)
# library(vegan)
# library(boral)

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
  

# #==================================================================
# # LATENT-VARIABLE (FACTOR ANALYSIS) MODELS FOR COVARIATES ALONE
# #==================================================================
# 
# #-------------------------------
# # JAGS
# #-------------------------------
# 
# # site-level covariates
# X <- as.matrix(lulc.roads.data[,c("ccap.hdev","ccap.open","ccap.lmdev","ccap.forest","ccap.wetland",
#                         "nlcd.imperv","restoration","pop.2010","roads.13","roads.4","roads.5",
#                         "roads.tot","traffic","area")])
# 
# # indices of covariates assumed to be beta (gamma) distributed
# beta.indx <- which(dimnames(X)[[2]] %in% c("ccap.open","ccap.lmdev","ccap.hdev","ccap.forest","ccap.wetland","nlcd.imperv"))
# gamma.indx <- which(dimnames(X)[[2]] %in% c("area","restoration","pop.2010","roads.13","roads.4","roads.5","roads.tot","traffic"))
# 
# # rescale gamma-distributed vars to SD = 1, bound all vars away from zero
# X[,gamma.indx] <- sweep(X[,gamma.indx], 2, apply(X[,gamma.indx], 2, sd), "/")
# X[X==0] <- 1e-4
# 
# S <- nrow(X)                     # number of sites
# D <- ncol(X)                     # dimension of site-level covariates
# L <- 1                           # dimension of latent-factor space (user-specified!)
# # n <- psm$n                       # number of females sampled
# # n.psm <- psm$n.psm               # number of PSM females
# # site <- as.numeric(psm$site)     # site index
# 
# # function to generate initial values for chains
# inits.fn <- function() list(a0 = rnorm(D, sapply(1:D, function(j) ifelse(j %in% beta.indx,
#                                                                          qlogis(mean(X[,j])),
#                                                                          log(mean(X[,j])))), 1),
#                             AA = matrix(rnorm(D*L, 0, 1), D, L),
#                             phi = sapply(1:D, function(j) ifelse(j %in% beta.indx, 
#                                                                    runif(1, 5, 10), 
#                                                                    runif(1, 0.5, 2))))
# 
# # Fit it!
# jags.lulc <- jags.parallel(data = c("beta.indx","gamma.indx","X","S","D","L"), 
#                          inits = inits.fn, 
#                          parameters.to.save = c("a0","A","Z","phi","mu"), 
#                          model.file = "beta_gamma_FA_JAGS.r",
#                          n.chains = 4, n.iter = 100000, n.burnin = 50000, n.thin = 50, 
#                          DIC = T, jags.seed = Sys.time())
# 
# # Plot and print summaries of results
# print(jags.lulc,2)
# 
# # Plot chains and histograms
# # plot(mcmc.list(mapply(function(d,x) mcmc(x$BUGS$sims.array[,d,]), 1:jags.psm$BUGS$n.chains, list(x=jags.psm), SIMPLIFY=F)), ask=T) 
# plot(mcmc.list(mapply(function(d,x) mcmc(x$BUGS$sims.array[,d,substring(dimnames(x$BUGS$sims.array)[[3]],1,2) != "mu"]), 
#                       1:jags.lulc$BUGS$n.chains, list(x=jags.lulc), SIMPLIFY=F)), ask=T)
# 
# # Fitted vs. observed
# par(mfrow=c(3,5))
# for(i in 1:D)
# {
#   plot(jags.lulc$BUGS$mean$mu[,i], X[,i], xlab="fitted", ylab="observed", main=dimnames(X)[[2]][i])
#   abline(0,1)
# }

# #-------------------------------
# # boral
# #-------------------------------
# 
# boral.lulc <- boral(y = X, family = ifelse(1:D %in% beta.indx, "beta", "gamma"), num.lv = 1,
#                     n.burnin = 5000, n.iteration = 35000, n.thin = 5,
#                     save.model = TRUE, seed = Sys.time(), calc.ics = TRUE,
#                     hypparams = c(100,100,100,1e5), ssvs.index = -1)
# 
# boral.lulc
# 
# # traceplots (not very helpful due to uninformative parameter names)
# plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1, list(x=boral.lulc$jags$BUGS$sims.array), SIMPLIFY=F)), ask=T)
# 
# # compare intercept, loading and factor score estimates from boral and my own JAGS code
# par(mfrow=c(1,3))
# plot(jags.lulc$BUGS$mean$a0, boral.lulc$jags$BUGS$mean$all.params[,1])
# plot(jags.lulc$BUGS$mean$A, boral.lulc$jags$BUGS$mean$all.params[,2])
# plot(jags.lulc$BUGS$mean$Z, boral.lulc$jags$BUGS$mean$lvs)
# 
# plot(jags.lulc$BUGS$mean$a0, boral.lulc$jags$BUGS$mean$all.params[,1]); abline(0,1)
# plot(jags.lulc$BUGS$mean$A, -boral.lulc$jags$BUGS$mean$all.params[,2]); abline(0,1)
# plot(jags.lulc$BUGS$mean$Z, -boral.lulc$jags$BUGS$mean$lvs); abline(0,1)


#==================================================================
# LATENT-VARIABLE REGRESSION MODELS FOR PSM
#==================================================================

#------------------------------------------------------
# Fit model to sites with PSM observations
#------------------------------------------------------

# site-level covariates

# composition data (multiple categories):
# bound away from 0 and 1 and re-standardize to sum to 1,
# then transform to log ratios
X1 <- as.matrix(lulc.roads.data[,c("ccap.hdev","ccap.ldev","ccap.mdev",
                                   "ccap.decid","ccap.evgrn","ccap.mxforest",
                                   "ccap.open","ccap.wetland","ccap.ag", "ccap.other")])
X1[X1==0] <- 1e-4
X1[X1==1] <- 1 - 1e-4
X1 <- sweep(X1, 1, rowSums(X1), "/")
X1 <- sweep(log(X1[,-ncol(X1)]), 1, log(X1[,ncol(X1)]), "-") 
X1 <- sweep(sweep(X1, 2, colMeans(X1), "-"), 2, apply(X1, 2, sd), "/")

# proportion data:
# bound away from 0 and 1 and logit-transform
X2 <- as.matrix(lulc.roads.data[,"nlcd.imperv",drop=F])
X2[X2==0] <- 1e-4
X2[X2==1] <- 1 - 1e-4
X2 <- qlogis(X2)
X2 <- sweep(sweep(X2, 2, colMeans(X2), "-"), 2, apply(X2, 2, sd), "/")

# nonnegative continuous data:
# scale to SD = 1, bound away from 0
X3 <- as.matrix(lulc.roads.data[,c("roads.1","roads.2","roads.3","roads.4","roads.5","traffic",
                                   "restoration","pop.census","pop.lscan")])
X3 <- sweep(X3, 2, apply(X3, 2, sd), "/")
X3[X3==0] <- 1e-4
X3 <- sweep(X3, 2, apply(X3, 2, sd), "/")

# reassemble all variables into matrix
X <- cbind(X2,X1,X3)

# indices of covariates assumed to be normally (gamma) distributed
normal.indx <- which(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd"))
gamma.indx <- which(!(substring(dimnames(X)[[2]],1,4) %in% c("ccap","nlcd")))

S <- nrow(X)                     # number of sites
D.normal <- length(normal.indx)  # dimension of normally distributed covariates
D.gamma <- length(gamma.indx)    # dimension of gamma-distributed covariates
D <- ncol(X)                     # dimension of site-level covariates
L <- 1                           # dimension of latent-factor space (user-specified!)
n <- psm$n                       # number of females sampled
n.psm <- psm$n.psm               # number of PSM females
site <- as.numeric(psm$site)     # site index

ppt.su <- as.vector(scale(psm$ppt.su/10, scale=F))    # center and convert mm to cm
ppt.fa <- as.vector(scale(psm$ppt.fa/10, scale=F))

# function to generate initial values for chains
inits.fn <- function() list(a0 = rnorm(D, c(colMeans(X[,1:D.normal]), colMeans(log(X[,-(1:D.normal)]))), 1),
                            AA = matrix(rnorm(D*L, 0, 1), D, L),
                            phi = runif(D, 0.5, 1),
                            B0 = rnorm(1, qlogis(mean(n.psm/n, na.rm=T)), 1),
                            B0.Z.px = rnorm(L, 0, 0.1),
                            sigma.b0 = runif(1, 0.1, 1),
                            B.su = rnorm(1, 0, 1),
                            B.su.Z.px = rnorm(L, 0, 0.1),
                            sigma.b.su = runif(1, 0.1, 1),
                            B.fa = rnorm(1, 0, 1),
                            B.fa.Z.px = rnorm(L, 0, 0.1),
                            sigma.b.fa = runif(1, 0.1, 1))

# Fit it!
system.time(
    jags.psm <- jags(data = c("normal.indx","gamma.indx","X","S","D","L","n","n.psm","site","ppt.su","ppt.fa"), 
                     inits = inits.fn, 
                     parameters.to.save = c("a0","A","Z","phi","mu",
                                            "B0","B0.Z","sigma.b0","b0",
                                            "B.su","B.su.Z","sigma.b.su","b.su",
                                            "B.fa","B.fa.Z","sigma.b.fa","b.fa",
                                            "p.psm"), 
                     model.file = "cohoPSM_SEM_JAGS.r",
                     n.chains = 3, n.iter = 100000, n.burnin = 50000, n.thin = 50, 
                     DIC = T, parallel = TRUE, seed = Sys.time())
)


# Plot and print summaries of results
print(jags.psm,2)

# Plot chains and histograms
plot(jags.psm)

# plot(mcmc.list(mapply(function(d,x) mcmc(x$BUGS$sims.array[,d,!substring(dimnames(x$BUGS$sims.array)[[3]],1,2) %in% c("mu","p.")]), 
#                       1:jags.psm$BUGS$n.chains, list(x=jags.psm), SIMPLIFY=F)), ask=T)


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
system.time(
  jags.psm.all <- jags.basic(data = c("beta.indx","dirch.indx","gamma.indx","X","S","D","L","n","n.psm","site","ppt.su","ppt.fa"), 
                           inits = inits.fn, 
                           parameters.to.save = c("a0","A","Z","phi","mu",
                                                  "B0","B0.Z","sigma.b0","b0",
                                                  "B.su","B.su.Z","sigma.b.su","b.su",
                                                  "B.fa","B.fa.Z","sigma.b.fa","b.fa",
                                                  "p.psm"), 
                           model.file = "cohoPSM_SEM_JAGS.r",
                           #                            n.chains = 3, n.iter = 100000, n.burnin = 20000, n.thin = 80, 
                           n.chains = 3, n.iter = 2000, n.burnin = 1000, n.thin = 1, 
                           DIC = T, parallel = TRUE, seed = Sys.time(), save.model = FALSE)
)

# Plot chains and histograms (if you dare)
plot(jags.psm.all)
plot(jags.psm.all[,substring(varnames(jags.psm.all),1,4)=="b.su"], ask=T)

# Store predicted P(PSM) in matrix
psm.pre <- as.matrix(jags.psm.all[,substring(varnames(jags.psm.all),1,5)=="p.psm"])
row.names(psm.pre) <- NULL
psm.pre <- data.frame(site=psm.all$site[psm.all$data=="pre"], 
                      psm.mean=colMeans(psm.pre)[psm.all$data=="pre"],
                      logit.psm.sd=apply(qlogis(psm.pre),2,sd)[psm.all$data=="pre"])
write.table(psm.pre, "PSM_predictions.txt", sep="\t", row.names=FALSE)


#-----------------------------------------
# Posteriors of loadings
#-----------------------------------------

dev.new(width=10, height=10)
# png(filename="loadings.png", width=10, height=10, units="in", res=300, type="cairo-png")
par(mar=c(5.1,9,2.1,2.1))

A <- jags.psm$sims.list$A[,-max(dirch.indx),1]  
bxp.dat <- boxplot(A, plot=F)
bxp.dat$stats[3,] <- colMeans(A)
bxp.dat$stats[c(2,4),] <- apply(A,2,quantile,c(0.05,0.95))
bxp.dat$stats[c(1,5),] <- apply(A,2,quantile,c(0.025,0.975))

bxp(bxp.dat, xlab="Factor loading (A)", ylab="", ylim=range(bxp.dat$stats), 
    yaxt="n", horizontal=T, boxwex=0.4, outpch="", whisklty=1, las=1, cex.axis=1.2, cex.lab=1.5)
axis(2, at=1:ncol(A), labels=dimnames(X)[[2]][-max(dirch.indx)], las=1, cex.axis=1.3)
abline(v=0, col="black", lty=2, lwd=2)
# dev.off()

#-----------------------------------------
# Observed PSM vs. factor scores
#-----------------------------------------

dev.new(width=10, height=10)
# png(filename="PSM_vs_Z.png", width=10, height=10, units="in", res=300, type="cairo-png")

plot(jags.psm$mean$Z[site], n.psm/n, xlab="Urbanization (Z)", ylab="Observed P(PSM)",
     pch="", cex.axis=1.2, cex.lab=1.5, las=1)
for(j in 1:S)
  curve(plogis(jags.psm$mean$b0[j] + jags.psm$mean$B0.Z*x), 
        from=-2, to=2, col=site[j], lwd=1, add=T)
points(jags.psm$mean$Z[site], n.psm/n, pch=16, col=site)
curve(plogis(jags.psm$mean$B0 + jags.psm$mean$B0.Z*x), from=-2, to=2, col="black", lwd=4, add=T)
# dev.off()

#----------------------------------------------------------------------------------------------
# Plots of site-level intercept and slope estimates against site-level LU/LC factor scores
#----------------------------------------------------------------------------------------------

dev.new(width=15,height=5)
# png(filename="psm_site-level_regressions.png", width=15, height=5, units="in", res=300, type="cairo-png")
par(mfrow=c(1,3), mar=c(5.1,5.1,4.1,2.1))

plot(jags.psm$mean$Z, jags.psm$mean$b0, pch="", las=1, cex.lab=1.8, cex.axis=1.5,
     xlim=range(apply(jags.psm$sims.list$Z, 2, quantile, c(0.025,0.975))),
     ylim=range(apply(jags.psm$sims.list$b0, 2, quantile, c(0.025,0.975))),
     xlab="Urbanization", ylab="Intercept")
fits <- sapply(1:jags.psm$mcmc.info$n.samples, 
               function(i,mod=jags.psm) 
                 mod$sims.list$B0[i] + mod$sims.list$B0.Z[i]*seq(par()$usr[1], par()$usr[2], length=50))
cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], alpha=0.5*255, maxColorValue=255)
polygon(c(seq(par()$usr[1], par()$usr[2], length=50), seq(par()$usr[2], par()$usr[1], length=50)),
        c(apply(fits,1,quantile,0.025), rev(apply(fits,1,quantile,0.975))),
        col=cc, border=NA)
abline(jags.psm$mean$B0, jags.psm$mean$B0.Z, lwd=2)
points(jags.psm$mean$Z, jags.psm$mean$b0, pch=18, las=1, 
       cex=1.5*log(table(site[!is.na(n.psm)])/min(table(site[!is.na(n.psm)])) + 1))
segments(jags.psm$mean$Z, apply(jags.psm$sims.list$b0, 2, quantile, 0.025), 
         jags.psm$mean$Z, apply(jags.psm$sims.list$b0, 2, quantile, 0.975))
segments(apply(jags.psm$sims.list$Z, 2, quantile, 0.025), jags.psm$mean$b0, 
         apply(jags.psm$sims.list$Z, 2, quantile, 0.975), jags.psm$mean$b0)
mtext("A", side=3, adj=0, cex=1.5)
legend("topleft", legend=c("1","2","5","10"), title="N years", y.intersp=1.5, pch=18, pt.cex=1.5*log(c(1,2,5,10) + 1))

plot(jags.psm$mean$Z, jags.psm$mean$b.su, pch="", las=1, cex.lab=1.8, cex.axis=1.5,
     xlim=range(apply(jags.psm$sims.list$Z, 2, quantile, c(0.025,0.975))),
     ylim=range(apply(jags.psm$sims.list$b.su, 2, quantile, c(0.025,0.975))),
     xlab="Urbanization", ylab=expression("Effect of summer precipitation (cm" * {}^-1 * ")"))
fits <- sapply(1:jags.psm$mcmc.info$n.samples, 
               function(i,mod=jags.psm) 
                 mod$sims.list$B.su[i] + mod$sims.list$B.su.Z[i]*seq(par()$usr[1], par()$usr[2], length=50))
cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], alpha=0.5*255, maxColorValue=255)
polygon(c(seq(par()$usr[1], par()$usr[2], length=50), seq(par()$usr[2], par()$usr[1], length=50)),
        c(apply(fits,1,quantile,0.025), rev(apply(fits,1,quantile,0.975))),
        col=cc, border=NA)
abline(jags.psm$mean$B.su, jags.psm$mean$B.su.Z, lwd=2)
points(jags.psm$mean$Z, jags.psm$mean$b.su, pch=18, las=1, 
       cex=1.5*log(table(site[!is.na(n.psm)])/min(table(site[!is.na(n.psm)])) + 1))
segments(jags.psm$mean$Z, apply(jags.psm$sims.list$b.su, 2, quantile, 0.025), 
         jags.psm$mean$Z, apply(jags.psm$sims.list$b.su, 2, quantile, 0.975))
segments(apply(jags.psm$sims.list$Z, 2, quantile, 0.025), jags.psm$mean$b.su, 
         apply(jags.psm$sims.list$Z, 2, quantile, 0.975), jags.psm$mean$b.su)
mtext("B", side=3, adj=0, cex=1.5)

plot(jags.psm$mean$Z, jags.psm$mean$b.fa, pch="", las=1, cex.lab=1.8, cex.axis=1.5,
     xlim=range(apply(jags.psm$sims.list$Z, 2, quantile, c(0.025,0.975))),
     ylim=range(apply(jags.psm$sims.list$b.fa, 2, quantile, c(0.025,0.975))),
     xlab="Urbanization", ylab=expression("Effect of fall precipitation (cm" * {}^-1 * ")"))
fits <- sapply(1:jags.psm$mcmc.info$n.samples, 
               function(i,mod=jags.psm) 
                 mod$sims.list$B.fa[i] + mod$sims.list$B.fa.Z[i]*seq(par()$usr[1], par()$usr[2], length=50))
cc <- col2rgb("darkgray")
cc <- rgb(cc[1], cc[2], cc[3], alpha=0.5*255, maxColorValue=255)
polygon(c(seq(par()$usr[1], par()$usr[2], length=50), seq(par()$usr[2], par()$usr[1], length=50)),
        c(apply(fits,1,quantile,0.025), rev(apply(fits,1,quantile,0.975))),
        col=cc, border=NA)
abline(jags.psm$mean$B.fa, jags.psm$mean$B.fa.Z, lwd=2)
points(jags.psm$mean$Z, jags.psm$mean$b.fa, pch=18, las=1, 
       cex=1.5*log(table(site[!is.na(n.psm)])/min(table(site[!is.na(n.psm)])) + 1))
segments(jags.psm$mean$Z, apply(jags.psm$sims.list$b.fa, 2, quantile, 0.025), 
         jags.psm$mean$Z, apply(jags.psm$sims.list$b.fa, 2, quantile, 0.975))
segments(apply(jags.psm$sims.list$Z, 2, quantile, 0.025), jags.psm$mean$b.fa, 
         apply(jags.psm$sims.list$Z, 2, quantile, 0.975), jags.psm$mean$b.fa)
mtext("C", side=3, adj=0, cex=1.5)
# dev.off()


#----------------------------------------------------------------------------
# Fit to data-rich, predict data-poor: plot predictions vs. observations 
#----------------------------------------------------------------------------

# dev.new()
png(filename="fit_data-rich_predict_data-poor.png", width=10, height=10, units="in", res=300, type="cairo-png")
plot(jags.psm$BUGS$mean$p.psm, psm$n.psm/psm$n, xlab="Predicted P(PSM)", ylab="Observed P(PSM)",
     pch="", cex.axis=1.2, cex.lab=1.5)
points(jags.psm$BUGS$mean$p.psm[!is.na(n.psm)], (psm$n.psm/psm$n)[!is.na(n.psm)], pch=0, cex=1.2, col="darkgray")
points(jags.psm$BUGS$mean$p.psm[is.na(n.psm)], (psm$n.psm/psm$n)[is.na(n.psm)], pch=16, cex=1.2, col="black")
abline(0,1)
legend("topleft", c("fitted","OOS prediction"), pch=c(0,16), cex=1.2, col=c("darkgray","black"))
# dev.off()




