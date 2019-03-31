




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

# add column of annual summer and fall precip to spawner data
ppt.su <- spatial.data[,names(spatial.data)=="site" | substring(names(spatial.data),1,6)=="ppt.su"]
ppt.su <- data.frame(site=rep(ppt.su$site, each=ncol(ppt.su) - 1), 
                     year=rep(substring(names(ppt.su[,-1]), 8), nrow(ppt.su)),
                     ppt.su=as.vector(t(as.matrix(ppt.su[,-1]))))
ppt.fa <- spatial.data[,names(spatial.data)=="site" | substring(names(spatial.data),1,6)=="ppt.fa"]
ppt.fa <- data.frame(site=rep(ppt.fa$site, each=ncol(ppt.fa) - 1), 
                     year=rep(substring(names(ppt.fa[,-1]), 8), nrow(ppt.fa)),
                     ppt.fa=as.vector(t(as.matrix(ppt.fa[,-1]))))
#----------------------------------------------


#----------------------------------------------
# ADDITIONAL SITES FOR PSM PREDICTIONS
#----------------------------------------------

# read in landscape data
spatial.data.pre <- read.csv("spatial.data.predict.csv", header=T)
spatial.data.pre <- spatial.data.pre[spatial.data.pre$area > 0.1 & rowSums(spatial.data.pre[,c(27:38,41)])!=0,]  # drop sites < 1 ha and "empty" sites

# reshape so each row is a site x year combination
ppt.su.pre <- spatial.data.pre[,names(spatial.data.pre)=="ID" | substring(names(spatial.data.pre),1,6)=="ppt.su"]
ppt.su.pre <- data.frame(ID=rep(ppt.su.pre$ID, each=ncol(ppt.su.pre) - 1), 
                         year=rep(substring(names(ppt.su.pre[,-1]), 8), nrow(ppt.su.pre)),
                         ppt.su=as.vector(t(as.matrix(ppt.su.pre[,-1]))))
ppt.fa.pre <- spatial.data.pre[,names(spatial.data.pre)=="ID" | substring(names(spatial.data.pre),1,6)=="ppt.fa"]
ppt.fa.pre <- data.frame(ID=rep(ppt.fa.pre$ID, each=ncol(ppt.fa.pre) - 1), 
                         year=rep(substring(names(ppt.fa.pre[,-1]), 8), nrow(ppt.fa.pre)),
                         ppt.fa=as.vector(t(as.matrix(ppt.fa.pre[,-1]))))

data.pre <- merge(ppt.su.pre, ppt.fa.pre, by=c("ID","year"))
#----------------------------------------------



#--------------------------------------------------------------
# Ordination for LU/LC and roads together
#--------------------------------------------------------------

# Data from PSM sampling sites
lulc.roads.data <- spatial.data[,c("ccap.open","ccap.lmdev","ccap.hdev","ccap.forest","ccap.wetland",
                                   "nlcd.imperv","restoration","pop.2010","roads.13","roads.4","roads.5",
                                   "roads.tot","traffic")]

# Data from unobserved sites
lulc.roads.data.pre <- spatial.data.pre[,c("ccap.open","ccap.lmdev","ccap.hdev","ccap.forest","ccap.wetland",
                                           "nlcd.imperv","restoration","pop.2010","roads.13","roads.4","roads.5",
                                           "roads.tot","traffic")]

# Combine PSM and unobserved data frames
lulc.roads.data.all <- rbind(data.frame(data="psm", lulc.roads.data), 
                             data.frame(data="pre", lulc.roads.data.pre))

# Transform variables
lulc.roads.data.all[,2:7] <- apply(lulc.roads.data.all[,2:7], 2, function(x) qlogis(pmin(pmax(x, 1e-5), 1 - 1e-6)))
lulc.roads.data.all[,8:14] <- apply(lulc.roads.data.all[,8:14], 2, function(x) log(x + min(x[x>0])))
# lulc.roads.data.all[,8:14] <- sweep(lulc.roads.data.all[,8:14], 2, apply(lulc.roads.data.all[,8:14], 2, max), "/")

# PCA
pca.lulc.roads.all <- prcomp(lulc.roads.data.all[lulc.roads.data.all$data=="psm",-1], scale=TRUE)
# pca.lulc.roads.all <- prcomp(lulc.roads.data.all[,-1], scale=TRUE)
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


# Merge spawner and landscape data frames
psm1 <- merge(spawner.data, ppt.su, by=c("site","year"))
psm1 <- merge(psm1, ppt.fa, by=c("site","year"))
psm1 <- merge(psm1, spatial.data[,substring(names(spatial.data),1,3) != "ppt"], by=c("ID","site"))
psm1 <- data.frame(psm1[,c("ID","site","watershed","source","year","n","n.psm","ppt.su","ppt.fa")], psm1[,-(1:9)])
psm1 <- psm1[order(psm1$site,psm1$year),]
psm <- psm1

# > names(psm)
# [1] "ID"             "site"           "watershed"      "source"         "year"          
# [6] "n"              "n.psm"          "ppt.su"         "ppt.fa"         "area"          
# [11] "ccap.open"      "ccap.lmdev"     "ccap.hdev"      "ccap.forest"    "ccap.wetland"  
# [16] "roads.13"       "roads.4"        "roads.5"        "roads.tot"      "nlcd.imperv"   
# [21] "restoration"    "pop.2010"       "pop.2000"       "pop.2011"       "traffic"       
# [26] "lulc.roads.pc1" "lulc.roads.pc2"

# Data for predictions
data.pre <- merge(data.pre, spatial.data.pre[,-(3:26)], by=c("ID"))
data.pre <- data.frame(data.pre, 
                       lulc.roads.pc1 = predict(pca.lulc.roads.all, 
                                                newdata=lulc.roads.data.all[lulc.roads.data.all$data=="pre",])[,1], 
                       lulc.roads.pc2 = predict(pca.lulc.roads.all, 
                                                newdata=lulc.roads.data.all[lulc.roads.data.all$data=="pre",])[,2])



#----------------------------------------------
# paired scatterplots of predictors
#----------------------------------------------
dev.new()
pairs(psm[,8:25], gap=0.2)
pairs(psm[,8:25][,!(names(psm[,8:25]) %in% c("ppt.su","ppt.fa","area","ccap.open","ccap.lmdev","roads.5","roads.tot","pop.2000","pop.2011"))])
pairs(psm[,c("ppt.su","ppt.fa","lulc.roads.pc1","lulc.roads.pc2")])


# strong correlations: 
# ccap.open vs. {ccap.lmdev, roads.tot, nlcd.imperv}
# ccap.lmdev vs. {ccap.forest, roads.4, roads.5, roads.tot, nldc.imperv, pop.2000, pop.2010, pop.2011}
# ccap.forest vs. {roads.4, roads.5, roads.tot}
# roads.4 vs. {roads.5, roads.tot, nlcd.imperv}
# roads.5 vs. {roads.tot, nlcd.imperv, pop.2000, pop.2010, pop.2011}
# roads.tot vs. {nlcd.imperv, pop.2000, pop.2010, pop.2011}
# nlcd.imperv vs. {pop.2000, pop.2010, pop.2011}
# pop.2000 vs. pop.2010 vs. pop.2011



#----------------------------------------------
# Assign variables to objects for JAGS
#----------------------------------------------
n <- psm$n
n.psm <- psm$n.psm
site <- as.numeric(psm$site)
X.dat <- data.frame(psm[,c("n.psm","n","ppt.su","ppt.fa","nlcd.imperv","roads.13","lulc.roads.pc1","lulc.roads.pc2")])
X.dat[,c("ppt.su","ppt.fa","nlcd.imperv","roads.13","lulc.roads.pc1","lulc.roads.pc2")] <- 
  scale(X.dat[,c("ppt.su","ppt.fa","nlcd.imperv","roads.13","lulc.roads.pc1","lulc.roads.pc2")])

# X <- model.matrix(n.psm/n ~ ppt.su + ppt.fa + lulc.roads.pc1 + lulc.roads.pc2 +
#                     ppt.su:lulc.roads.pc1 + ppt.su:lulc.roads.pc2 +
#                     ppt.fa:lulc.roads.pc1 + ppt.fa:lulc.roads.pc2,
#                   data = X.dat)

# X <- model.matrix(n.psm/n ~ ppt.su + ppt.fa + nlcd.imperv + roads.13 +
#                     ppt.su:nlcd.imperv + ppt.su:roads.13 +
#                     ppt.fa:nlcd.imperv + ppt.fa:roads.13,
#                   data = X.dat)

X <- model.matrix(n.psm/n ~ poly(ppt.su,2,raw=T) + poly(ppt.fa,2,raw=T) + lulc.roads.pc1 +
                    ppt.su:ppt.fa + ppt.su:lulc.roads.pc1 + ppt.fa:lulc.roads.pc1,
                  data = X.dat)

X.dat.pre <- data.frame(year=2011, data.pre[data.pre$year==2011,c("ID","ppt.su","ppt.fa","nlcd.imperv","roads.13","lulc.roads.pc1","lulc.roads.pc2")])
X.dat.pre[,-(1:2)] <- as.data.frame(scale(X.dat.pre[,-(1:2)], 
                                          center=colMeans(psm[,names(X.dat.pre)[-(1:2)]]), 
                                          scale=apply(psm[,names(X.dat.pre)[-(1:2)]], 2, sd)))

# X.pre <- model.matrix(rep(1,nrow(X.dat.pre)) ~ ppt.su + ppt.fa + lulc.roads.pc1 + lulc.roads.pc2 +
#                         ppt.su:lulc.roads.pc1 + ppt.su:lulc.roads.pc2 +
#                         ppt.fa:lulc.roads.pc1 + ppt.fa:lulc.roads.pc2,
#                   data = X.dat.pre)

# X.pre <- model.matrix(rep(1,nrow(X.dat.pre)) ~ ppt.su + ppt.fa + nlcd.imperv + roads.13 +
#                         ppt.su:nlcd.imperv + ppt.su:roads.13 +
#                         ppt.fa:nlcd.imperv + ppt.fa:roads.13,
#                       data = X.dat.pre)

X.pre <- model.matrix(rep(1,nrow(X.dat.pre)) ~ poly(ppt.su,2,raw=T) + poly(ppt.fa,2,raw=T) + lulc.roads.pc1 +
                        ppt.su:ppt.fa + ppt.su:lulc.roads.pc1 + ppt.fa:lulc.roads.pc1,
                      data = X.dat.pre)



#==================================================================
# MODELS OF PRESPAWNING MORTALITY
#==================================================================

#----------------------------------------------
# GLM, no random effects
#----------------------------------------------

glm.psm <- glm((n.psm/n) ~ ppt.su + ppt.fa + nlcd.imperv + roads.13 +
                 ppt.su:nlcd.imperv + ppt.su:roads.13 +
                 ppt.fa:nlcd.imperv + ppt.fa:roads.13,
               data = X.dat, family = binomial(link="logit"), weights = n)

# glm.psm <- glm((n.psm/n) ~ ppt.su + ppt.fa + lulc.roads.pc1 + lulc.roads.pc2 +
#                  ppt.su:lulc.roads.pc1 + ppt.su:lulc.roads.pc2 +
#                  ppt.fa:lulc.roads.pc1 + ppt.fa:lulc.roads.pc2,
#                data = X.dat, family = binomial(link="logit"), weights = n)

summary(glm.psm)


#----------------------------------------------
# GLMM via lme4, random site-level intercept
#----------------------------------------------

# glmer.psm <- glmer((n.psm/n) ~  (1|site) + ppt.su + ppt.fa + lulc.roads.pc1 + lulc.roads.pc2 +
#                      ppt.su:lulc.roads.pc1 + ppt.su:lulc.roads.pc2 +
#                      ppt.fa:lulc.roads.pc1 + ppt.fa:lulc.roads.pc2,
#                    data = X.dat, family = binomial(link="logit"), weights = n, 
#                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e5)))

# glmer.psm <- glmer((n.psm/n) ~  (1|site) + ppt.su + ppt.fa + nlcd.imperv + roads.13 +
#                      ppt.su:nlcd.imperv + ppt.su:roads.13 +
#                      ppt.fa:nlcd.imperv + ppt.fa:roads.13,
#                    data = X.dat, family = binomial(link="logit"), weights = n, 
#                    glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e5)))

glmer.psm <- glmer((n.psm/n) ~  (1|site) + poly(ppt.su,2,raw=T) + poly(ppt.fa,2,raw=T) + lulc.roads.pc1 +
                     ppt.su:ppt.fa + ppt.su:lulc.roads.pc1 + ppt.fa:lulc.roads.pc1,
                   data = X.dat, family = binomial(link="logit"), weights = n, 
                   glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e5)))

summary(glmer.psm)


#----------------------------------------------
# Bayesian GLMM via JAGS
#----------------------------------------------

inits.fn <- function() list(b = rnorm(ncol(X), fixef(glmer.psm), 1), sigma.b0.site = runif(1,0.5,2))

jags.psm <- jags(data = c("X", "site", "n", "n.psm", "X.pre"), 
                 inits = inits.fn, 
                 #                  parameters.to.save = c("b", "b0.site", "sigma.b0.site","lp.psm"), 
                 parameters.to.save = c("b", "b0.site", "sigma.b0.site","lp.psm", "lp.psm.pre"), 
                 model.file = "PSM_GLMM_JAGS.r",
                 n.chains = 3, n.iter = 100000, n.burnin = 5000, n.thin = 95, 
                 DIC = T, jags.seed = Sys.time())

# Plot and print summaries of results
print(jags.psm,2)

# Plot chains and histograms
plot(mcmc.list(mapply(function(d,x) mcmc(x[,d,]), 1:3, list(x=jags.psm$BUGS$sims.array), SIMPLIFY=F)), ask=T)

# QQ plot of random effects
dev.new()
qqnorm(jags.psm$BUGS$mean$b0.site)
qqline(jags.psm$BUGS$mean$b0.site)

# Paired scatterplots of posteriors for fixed effects
dev.new(width=10,height=10)
pairs(jags.psm$BUGS$sims.list$b, labels=dimnames(X)[[2]])

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


# Add predicted P(PSM) to matrix of predicted sites
X.dat.pre <- data.frame(X.dat.pre, p.psm = colMeans(plogis(jags.psm$BUGS$sims.list$lp.psm.pre)),
                        sd.logit.psm = apply(jags.psm$BUGS$sims.list$lp.psm.pre, 2, sd),
                        sd.p.psm = apply(plogis(jags.psm$BUGS$sims.list$lp.psm.pre), 2, sd))

write.table(X.dat.pre, "PSM_predictions.txt", row.names=F, sep="\t")



#==================================================================
# ANCILLARY PLOTS
#==================================================================

# Pairwise scatterplots of all predictors (main effects and intxns) and predicted P(PSM)
dev.new(width=14,height=14)
pairs(data.frame(X.pre[,-1], lp.psm=jags.psm$BUGS$mean$lp.psm.pre))






