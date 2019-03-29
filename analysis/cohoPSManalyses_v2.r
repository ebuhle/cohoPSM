




options(device=windows)
library(lme4)
# library(MCMCglmm)
library(R2jags)
library(vegan)


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


#----------------------------------------------
# ADDITIONAL SITES FOR PSM PREDICTIONS
#----------------------------------------------

# read in landscape data
spatial.data.pre <- read.csv("spatial.data.predict.csv", header=T)
spatial.data.pre <- spatial.data.pre[spatial.data.pre$area > 0.1 &  # drop sites < 1 ha and "empty" sites
                                       rowSums(spatial.data.pre[,!(substring(names(spatial.data.pre),1,2) %in% c("ID","si","wa","ar","pp"))])!=0,] 


#--------------------------------------------------------------
# Ordination for LU/LC and roads together
#--------------------------------------------------------------

# Data from PSM sampling sites
lulc.roads.data <- spatial.data[,c("ccap.open","ccap.lmdev","ccap.hdev","ccap.forest","ccap.wetland",
                             "nlcd.imperv","restoration","pop.2010","roads.13","roads.4","roads.5",
                             "roads.tot","traffic")]

# Data from unobserved sites
lulc.roads.data.pre <- spatial.data.pre[,names(lulc.roads.data)]

# Combine PSM and unobserved data frames
lulc.roads.data.all <- rbind(data.frame(data="psm", lulc.roads.data), 
                             data.frame(data="pre", lulc.roads.data.pre))

# Transform variables
lulc.roads.data.all[,2:7] <- apply(lulc.roads.data.all[,2:7], 2, function(x) qlogis(pmin(pmax(x, 1e-5), 1 - 1e-6)))
lulc.roads.data.all[,8:14] <- apply(lulc.roads.data.all[,8:14], 2, function(x) log(x + min(x[x>0])))
# lulc.roads.data.all[,8:14] <- sweep(lulc.roads.data.all[,8:14], 2, apply(lulc.roads.data.all[,8:14], 2, max), "/")

# PCA
# pca.lulc.roads.all <- prcomp(lulc.roads.data.all[lulc.roads.data.all$data=="psm",-1], scale=TRUE)
pca.lulc.roads.all <- prcomp(lulc.roads.data.all[,-1], scale=TRUE)
pca.lulc.roads.all$rotation[,1:2] <- -pca.lulc.roads.all$rotation[,1:2]  # flip sign so PC1, PC2 are pos. corr. w/ development
pca.lulc.roads.all$x[,1:2] <- -pca.lulc.roads.all$x[,1:2] 
pca.lulc.roads.all
summary(pca.lulc.roads.all)

dev.new(width=14,height=10)  # all pairwise biplots
par(mfrow=c(2,3), mar=c(3.1,4.1,1.1,1.1))
for(i in 1:3)
  for(j in (i+1):4)
    biplot(pca.lulc.roads.all, choices=c(i,j), col=c("gray","red"), cex=1.2)

dev.new()
# png(filename="LULC_roads_PCA.png", width=7, height=7, units="in", res=200, type="cairo-png")
biplot(pca.lulc.roads.all, las=1, cex.lab=1.5, col=c("gray","red"), cex.axis=1.2, cex=0.8)
# dev.off()

# # NMDS
# mds.lulc.roads.all <- metaMDS(lulc.roads.data.all[,-1], k=2, 
#                               trymax=20, trace=2, autotransform=F, wascores=T, noshare=F,
#                               maxit=5000, sratmax=1 - 1e-6, plot=T)
# mds.lulc.roads.all <- metaMDS(lulc.roads.data.all[,-1], k=2,
#                               previous.best = mds.lulc.roads.all,
#                               trymax=20, trace=2, autotransform=F, wascores=T, noshare=F,
#                               maxit=5000, sratmax=1 - 1e-6, plot=T)
# mds.lulc.roads.all$points <- -mds.lulc.roads.all$points     # reverse so axes go from less to more developed
# mds.lulc.roads.all$species <- -mds.lulc.roads.all$species
# mds.lulc.roads.all
# 
# dev.new(width=14,height=10)
# par(mfrow=c(1,2))
# stressplot(mds.lulc.roads.all, las=1)
# plot(mds.lulc.roads.all, las=1, type="n")
# points(mds.lulc.roads.all, display="sites", choices=c(1,2), pch=1, col="darkgray")
# text(mds.lulc.roads.all, display="species", choices=c(1,2), col="blue", cex=0.8, 
#      labels=c(ifelse(substring(names(lulc.roads.data.all[-1]), 1, 5) %in% c("ccap.","nlcd."),
#                      substring(names(lulc.roads.data.all[-1]), 6), names(lulc.roads.data.all[-1]))))
# legend("topright", legend=paste("stress =", round(mds.lulc.roads.all$stress, 2)), bty="n")



#--------------------------------------------------------------
# Add LU/LC + roads combined PC1 and PC2 to spatial data
#--------------------------------------------------------------

# Spatial data for PSM observations
spatial.data <- data.frame(spatial.data, 
                           lulc.roads.pc1 = predict(pca.lulc.roads.all, 
                                                    newdata=lulc.roads.data.all[lulc.roads.data.all$data=="psm",])[,1], 
                           lulc.roads.pc2 = predict(pca.lulc.roads.all, 
                                                    newdata=lulc.roads.data.all[lulc.roads.data.all$data=="psm",])[,2])
# Data for predictions
spatial.data.pre <- data.frame(spatial.data.pre, 
                       lulc.roads.pc1 = predict(pca.lulc.roads.all, 
                                                newdata=lulc.roads.data.all[lulc.roads.data.all$data=="pre",])[,1], 
                       lulc.roads.pc2 = predict(pca.lulc.roads.all, 
                                                newdata=lulc.roads.data.all[lulc.roads.data.all$data=="pre",])[,2])


#--------------------------------------------------------------
# Merge spawner and landscape data frames for PSM samples
#--------------------------------------------------------------

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
spatial.data.columns <- substring(names(spatial.data),1,3) != "ppt"
psm <- merge(psm, spatial.data[,spatial.data.columns], by=c("ID","site"))
psm <- data.frame(psm[,c("ID","site","watershed","source","year","n","n.psm","ppt.su","ppt.fa","area")], 
                  psm[,-(1:10)])
psm <- psm[order(psm$site,psm$year),]


#----------------------------------------------
# paired scatterplots of predictors
#----------------------------------------------

# PSM data
dev.new()
pairs(psm[,names(lulc.roads.data.all)[-1]], gap=0.2)

dev.new()
pairs(psm[,c("ppt.su","ppt.fa","lulc.roads.pc1","lulc.roads.pc2")])

# predictive data
dev.new()
pairs(spatial.data.pre[,names(lulc.roads.data.all)[-1]], gap=0.2)

dev.new()
pairs(spatial.data.pre[,c("ppt.su.2011","ppt.fa.2011","lulc.roads.pc1","lulc.roads.pc2")])


#----------------------------------------------
# Assign variables to objects for JAGS
#----------------------------------------------

#### PSM DATA
n <- psm$n
n.psm <- psm$n.psm
site <- as.numeric(psm$site)

# site-level variables
lulc.vars <- spatial.data[,c("area","ccap.open","ccap.lmdev","ccap.hdev","ccap.forest","ccap.wetland",
                                   "nlcd.imperv","restoration","pop.2010","roads.13","roads.4","roads.5",
                                   "roads.tot","traffic")]
# lulc.vars[,1:6] <- apply(lulc.vars[,1:6], 2, function(x) qlogis(pmin(pmax(x, 1e-5), 1 - 1e-6)))
# lulc.vars[,7:13] <- apply(lulc.vars[,7:13], 2, function(x) log(x + min(x[x>0])))
lulc.vars <- as.matrix(scale(lulc.vars))
lulc.ref <- prcomp(lulc.vars, scale = T)$rotation[,"PC1"]
lulc.ref <- lulc.ref*sign(lulc.ref["nlcd.imperv"])
# lulc <- as.vector(scale(spatial.data$lulc.roads.pc1))

# observation-level variables
ppt.su <- as.vector(scale(psm$ppt.su/10, scale=F))    # center and convert mm to cm
ppt.fa <- as.vector(scale(psm$ppt.fa/10, scale=F))


#### DATA FOR PREDICTIONS
year.pre <- 2011   # pick a year

# site-level variables
lulc.pre <- (spatial.data.pre$lulc.roads.pc1 - mean(psm$lulc.roads.pc1))/sd(psm$lulc.roads.pc1)

# observation-level variables
# (set precip anomalies to zero if you want predictions for an "average" year)
ppt.su.pre <- as.vector(scale(spatial.data.pre[,paste("ppt.su", year.pre, sep=".")]/10, 
                              center=mean(psm$ppt.su), scale=F))  # center (to mean in PSM dataset) and convert mm to cm
ppt.fa.pre <- as.vector(scale(spatial.data.pre[,paste("ppt.fa", year.pre, sep=".")]/10, 
                              center=mean(psm$ppt.fa), scale=F))  # center (to mean in PSM dataset) and convert mm to cm


#==================================================================
# MODELS OF PRESPAWNING MORTALITY
#==================================================================

#----------------------------------------------
# GLMM via lme4, random site-level intercept
#----------------------------------------------

glmer.dat <- data.frame(site = psm$site, 
                        ppt.su = scale(psm$ppt.su/10, scale=F), 
                        ppt.fa = scale(psm$ppt.fa/10, scale=F), 
                        lulc=scale(psm$lulc.roads.pc1), lulc.roads.pc2=scale(psm$lulc.roads.pc2))

glmer.psm <- glmer((n.psm/n) ~  (1|site) + ppt.su + ppt.fa + lulc +
                     ppt.su:lulc + ppt.fa:lulc,
                   data = glmer.dat, family = binomial(link="logit"), weights = n, 
                   glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e5)))

summary(glmer.psm)


#----------------------------------------------
# Bayesian GLMM via JAGS
#----------------------------------------------

inits.fn <- function() list(C = rnorm(ncol(lulc.vars),0,0.1),
                            b0.0 = rnorm(1, qlogis(mean(n.psm/n)), 1), b0.lulc = 0, sigma.b0 = runif(1,0.5,2),
                            b0 = qlogis(pmin(pmax(tapply(n.psm/n, site, mean), 0.001), 0.999)),
                            b.ppt.su.0 = rnorm(1,0,0.1), b.ppt.su.lulc = rnorm(1,0,0.1), sigma.b.ppt.su = runif(1,0.5,2),
                            b.ppt.su = rnorm(length(unique(site)),0,0.1),
                            b.ppt.fa.0 = rnorm(1,0,0.1), b.ppt.fa.lulc = rnorm(1,0,0.1), sigma.b.ppt.fa = runif(1,0.5,2),
                            b.ppt.fa = rnorm(length(unique(site)),0,0.1))
#                             b.ppt.su.fa.0 = rnorm(1,0,0.1), b.ppt.su.fa.lulc = 0, sigma.b.ppt.su.fa = runif(1,0.5,2))

jags.psm <- jags.parallel(data = c("ppt.su", "ppt.fa", "lulc.vars", "lulc.ref", "site", "n", "n.psm", "ppt.su.pre", "ppt.fa.pre", "lulc.pre"), 
                          inits = inits.fn, 
                          parameters.to.save = c("c", "lulc", "b0.0", "b0.lulc", "b0.lulc.effect", "sigma.b0", "b0",
                                                 "b.ppt.su.0", "b.ppt.su.lulc", "b.ppt.su.lulc.effect", "sigma.b.ppt.su", "b.ppt.su",
                                                 "b.ppt.fa.0", "b.ppt.fa.lulc", "b.ppt.fa.lulc.effect", "sigma.b.ppt.fa", "b.ppt.fa"),
                          #                                         "b.ppt.su.fa.0", "b.ppt.su.fa.lulc", "sigma.b.ppt.su.fa", "b.ppt.su.fa"),
                          #                                         "lp.psm"), 
                          #                                         "lp.psm", "lp.psm.pre"), 
                          model.file = "PSM_GLMM_JAGS.r",
                          n.chains = 3, n.iter = 500000, n.burnin = 100000, n.thin = 200, 
                          DIC = T, jags.seed = Sys.time())


# Add predicted P(PSM) to matrix of predicted sites
X.dat.pre <- data.frame(ID=spatial.data.pre$ID, 
                        p.psm = colMeans(plogis(jags.psm$BUGS$sims.list$lp.psm.pre)),
                        sd.logit.psm = apply(jags.psm$BUGS$sims.list$lp.psm.pre, 2, sd),
                        sd.p.psm = apply(plogis(jags.psm$BUGS$sims.list$lp.psm.pre), 2, sd))

write.table(X.dat.pre, "PSM_predictions.txt", row.names=F, sep="\t")

dev.new()
pairs(data.frame(ppt.su.pre=ppt.su.pre, ppt.fa.pre=ppt.fa.pre, lulc=lulc.pre, lp.psm.pre=jags.psm$BUGS$mean$lp.psm.pre))


# Plot and print summaries of results
print(jags.psm,2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1:3, list(x=jags.psm$BUGS$sims.array), SIMPLIFY=F)), ask=T)

# QQ plot of site-level estimates
dev.new(width=15, height=5)
par(mfrow=c(1,3))
qqnorm(jags.psm$BUGS$mean$b0); qqline(jags.psm$BUGS$mean$b0)             # site-level intercept
qqnorm(jags.psm$BUGS$mean$b.ppt.su); qqline(jags.psm$BUGS$mean$b.ppt.su) # site-level summer ppt slope
qqnorm(jags.psm$BUGS$mean$b.ppt.fa); qqline(jags.psm$BUGS$mean$b.ppt.fa) # site-level fall ppt slope

# Fitted vs. observed P(PSM), symbols indicate sites
repeats <- as.numeric(names(which(table(site) > 1)))  # 22 sites with >1 obs
rep.pch <- rep(2:12,2)
rep.col <- rep(1:11, each=2)

dev.new()
plot(colMeans(plogis(jags.psm$BUGS$sims.list$lp.psm)), psm$n.psm/psm$n, xlim=c(0,1), ylim=c(0,1),
     xlab="Fitted (posterior mean)", ylab="Observed", main="Probability of PSM (symbols = sites)", 
     cex.axis=1.2, cex.lab=1.5, las=1, pch="")
for(i in 1:length(levels(psm$site)))
 points(colMeans(plogis(jags.psm$BUGS$sims.list$lp.psm))[site==i], (psm$n.psm/psm$n)[site==i],
        pch=ifelse(i %in% repeats, rep.pch[which(repeats==i)], 1), 
        col=ifelse(i %in% repeats, rep.col[which(repeats==i)], 1))
abline(0,1)

# Fitted vs. observed P(PSM), panel for each year
dev.new(width=12, height=10)
par(mfrow=c(3,4))
for(i in sort(unique(psm$year)))
{
plot(colMeans(plogis(jags.psm$BUGS$sims.list$lp.psm))[psm$year==i], (psm$n.psm/psm$n)[psm$year==i], 
     xlim=c(0,1), ylim=c(0,1), cex.axis=1.2, cex.lab=1.5, las=1,
     xlab="Fitted (posterior mean)", ylab="Observed", main=paste("Probability of PSM (", i, ")", sep=""))
abline(0,1)
}

sapply(levels(factor(psm$year)), function(i) cor((psm$n.psm/psm$n)[psm$year==i], 
                                         colMeans(plogis(jags.psm$BUGS$sims.list$lp.psm))[psm$year==i]))


# Paired scatterplots of posteriors for "fixed" effects
dev.new(width=10,height=10)
pairs(data.frame(b0.0=jags.psm$BUGS$sims.list$b0.0, b0.lulc=jags.psm$BUGS$sims.list$b0.lulc,
                 b.ppt.su.0=jags.psm$BUGS$sims.list$b.ppt.su.0, b.ppt.su.lulc=jags.psm$BUGS$sims.list$b.ppt.su.lulc,
                 b.ppt.fa.0=jags.psm$BUGS$sims.list$b.ppt.fa.0, b.ppt.fa.lulc=jags.psm$BUGS$sims.list$b.ppt.fa.lulc))

# Plots of site-level intercept and slope estimates against site-level LU/LC covariate
dev.new(width=15,height=5)
# png(filename="psm_site-level_regressions.png", width=15, height=5, units="in", res=300, type="cairo-png")
par(mfrow=c(1,3), mar=c(5.1,5.1,4.1,2.1))

plot(lulc, jags.psm$BUGS$mean$b0, pch=16, las=1, cex=1.2, cex.lab=1.8, cex.axis=1.5,
     ylim=range(apply(jags.psm$BUGS$sims.list$b0, 2, quantile, c(0.025,0.975))),
     xlab="Land use/land cover PC1", ylab="Intercept")
segments(lulc, apply(jags.psm$BUGS$sims.list$b0, 2, quantile, 0.025), 
         lulc, apply(jags.psm$BUGS$sims.list$b0, 2, quantile, 0.975))
abline(jags.psm$BUGS$mean$b0.0, jags.psm$BUGS$mean$b0.lulc)
mtext("A", side=3, adj=0, cex=1.5)

plot(lulc, jags.psm$BUGS$mean$b.ppt.su, pch=16, las=1, cex=1.2, cex.lab=1.8, cex.axis=1.5,
     ylim=range(apply(jags.psm$BUGS$sims.list$b.ppt.su, 2, quantile, c(0.025,0.975))),
     xlab="Land use/land cover PC1", ylab=expression("Effect of summer precipitation (cm" * {}^-1 * ")"))
segments(lulc, apply(jags.psm$BUGS$sims.list$b.ppt.su, 2, quantile, 0.025), 
         lulc, apply(jags.psm$BUGS$sims.list$b.ppt.su, 2, quantile, 0.975))
abline(jags.psm$BUGS$mean$b.ppt.su.0, jags.psm$BUGS$mean$b.ppt.su.lulc)
mtext("B", side=3, adj=0, cex=1.5)

plot(lulc, jags.psm$BUGS$mean$b.ppt.fa, pch=16, las=1, cex=1.2, cex.lab=1.8, cex.axis=1.5,
     ylim=range(apply(jags.psm$BUGS$sims.list$b.ppt.fa, 2, quantile, c(0.025,0.975))),
     xlab="Land use/land cover PC1", ylab=expression("Effect of fall precipitation (cm" * {}^-1 * ")"))
segments(lulc, apply(jags.psm$BUGS$sims.list$b.ppt.fa, 2, quantile, 0.025), 
         lulc, apply(jags.psm$BUGS$sims.list$b.ppt.fa, 2, quantile, 0.975))
abline(jags.psm$BUGS$mean$b.ppt.fa.0, jags.psm$BUGS$mean$b.ppt.fa.lulc)
mtext("C", side=3, adj=0, cex=1.5)

# dev.off()



# Boxplot of marginal posteriors for fixed effects

# X.prettynames <- c(expression(intercept), expression("summer ppt"), expression("fall ppt"), 
#                    expression(PC1), expression(PC2), 
#                    expression("summer ppt " %*% " PC1"), expression("summer ppt " %*% " PC2"),
#                    expression("fall ppt " %*% " PC1"), expression("fall ppt " %*% " PC2"))
X.prettynames <- c(expression(intercept), expression("summer ppt"), expression("fall ppt"), 
                   expression("proportion impervious"), expression("class 1-3 road density"), 
                   expression("summer ppt " %*% " imperv"), expression("summer ppt " %*% " roads"),
                   expression("fall ppt " %*% " imperv"), expression("fall ppt " %*% " roads"))
X.cols <- c("royalblue","sienna","black","darkgray","royalblue","royalblue","sienna","sienna")
X.fill <- c("royalblue","sienna","black","darkgray",0,0,0,0)
X.medcol <- c("white","white","white","white","black","darkgray","black","darkgray")

bmat <- jags.psm$BUGS$sims.list$b        
bmat <- bmat[,-1]; labs <- X.prettynames[-1]   # don't show intercept

bxp.dat <- boxplot(bmat, plot=F)
bxp.dat$stats[3,] <- colMeans(bmat)
bxp.dat$stats[c(2,4),] <- apply(bmat,2,quantile,c(0.05,0.95))
bxp.dat$stats[c(1,5),] <- apply(bmat,2,quantile,c(0.025,0.975))
bxp.dat$stats <- bxp.dat$stats[,ncol(bxp.dat$stats):1]  # reverse order

dev.new()
# png(filename="psm_coefs.png", width=7, height=7, units="in", res=300, type="cairo-png")
par(mar=c(5.1,9,5.1,2.1))

bxp(bxp.dat, xlab=expression(paste("Standardized coefficient (", beta, ")")), ylab="", ylim=range(bxp.dat$stats), 
    yaxt="n", horizontal=T, whisklty=0, staplewex=0, boxwex=0.2, outpch="", 
    boxcol=0, boxfill=0, whiskcol=0, medlty="blank", medpch="", medcex=1, 
    las=1, cex.axis=1, cex.lab=1.2)
axis(2, at=1:ncol(bmat), labels=rev(labs), las=1)
abline(v=0, col="black", lty=2, lwd=2)
bxp(bxp.dat, horizontal=T, whisklty=1, staplewex=0, boxwex=0.2, outpch="", xlab="", ylab="", xaxt="n", yaxt="n",
    boxcol=rev(X.cols), boxfill=rev(X.fill), whiskcol=rev(X.cols), medlty="blank", medpch=16, medcex=1.2, 
    medcol=rev(X.medcol), medbg=rev(X.medcol), add=T)
arrows(c(-0.05, 0.05), rep(par()$usr[4] + 0.5, 2), c(par()$usr[1:2]), rep(par()$usr[4] + 0.5, 2), xpd=T, lwd=5)
mtext(c("decreasing PSM","increasing PSM"), side=3, line=2.5, at=c(-0.5, 1), adj=0.5)
rm(list=c("X.prettynames","X.cols","X.fill","X.medcol","bxp.dat","bmat","labs"))
# dev.off()



#==================================================================
# ANCILLARY PLOTS
#==================================================================

# Pairwise scatterplots of all predictors (main effects and intxns) and predicted P(PSM)
dev.new(width=14,height=14)
pairs(data.frame(X.pre[,-1], lp.psm=jags.psm$BUGS$mean$lp.psm.pre))

# Compare temporal and spatial variability in summer and fall precip in PSM DATA
ppt.data2 <- spatial.data[,names(spatial.data) %in% c("site",paste("ppt.su",2000:2011,sep="."),paste("ppt.fa",2000:2011,sep="."))]
ppt.su.cols <- names(ppt.data2) %in% paste("ppt.su",2000:2011,sep=".")
ppt.fa.cols <- names(ppt.data2) %in% paste("ppt.fa",2000:2011,sep=".")
ppt.data2 <- data.frame(site=rep(ppt.data2$site, each=n.yrs.ppt), 
                        year=rep(2000:2011, nrow(ppt.data2)),
                        ppt.su=as.vector(t(as.matrix(ppt.data2[,ppt.su.cols]))),
                        ppt.fa=as.vector(t(as.matrix(ppt.data2[,ppt.fa.cols]))))

lm.ppt.su <- lm(log(ppt.su) ~ site + factor(year), data=ppt.data2)
summary(lm.ppt.su)
anova(lm.ppt.su)

lm.ppt.fa <- lm(log(ppt.fa) ~ site + factor(year), data=ppt.data2)
summary(lm.ppt.fa)
anova(lm.ppt.fa)

dev.new()
par(mfrow=c(2,1))
plot(ppt.data2$year, ppt.data2$ppt.su, pch="", xlab="year", ylab="summer precip")
for(i in levels(ppt.data2$site))
  lines(ppt.data2$year[ppt.data2$site==i], ppt.data2$ppt.su[ppt.data2$site==i])

plot(ppt.data2$year, ppt.data2$ppt.fa, pch="", xlab="year", ylab="fall precip")
for(i in levels(ppt.data2$site))
  lines(ppt.data2$year[ppt.data2$site==i], ppt.data2$ppt.fa[ppt.data2$site==i])

# Compare temporal and spatial variability in summer and fall precip in PREDICTIVE DATA
dev.new()
par(mfrow=c(2,1))
cc <- col2rgb("darkgray")
cc <- rgb(red=cc[1], green=cc[2], blue=cc[3], alpha=255*0.5, maxColorValue=255)
plot(rep(2000:2011, each=nrow(spatial.data.pre)), 
     as.vector(as.matrix(spatial.data.pre[,paste("ppt.su", 2000:2011, sep=".")])),
     pch="", xlab="year", ylab="summer precip")
for(i in 1:nrow(spatial.data.pre))
  lines(2000:2011, unlist(spatial.data.pre[i,paste("ppt.su", 2000:2011, sep=".")]), col=cc)

plot(rep(2000:2011, each=nrow(spatial.data.pre)), 
     as.vector(as.matrix(spatial.data.pre[,paste("ppt.fa", 2000:2011, sep=".")])),
     pch="", xlab="year", ylab="fall precip")
for(i in 1:nrow(spatial.data.pre))
  lines(2000:2011, unlist(spatial.data.pre[i,paste("ppt.fa", 2000:2011, sep=".")]), col=cc)





