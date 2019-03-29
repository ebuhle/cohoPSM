#======================================================
# COHO PRESPAWN MORTALITY ANALYSIS
# (Feist et al. 2010, Science, in press)
#
# Ordination of catchment LU/LC data, and mixed
# models of mortality using raw layer data as
# predictors (ordination = no good).
#======================================================

#------------------------------------------------------
# DATA IMPORT AND MANIPULATION
#------------------------------------------------------

# Import main catchment dataset 
# 6 study sites are included, overlapping catchments are excluded,
# and all catchments < 1e6 ft2 are excluded
catchments.final <- read.table("catchments.final.txt","\t",header=T)

# Convert Impervious to a "proportion" scale
catchments.final$Impervious <- catchments.final$Impervious/100

# Import precipitation data
precip <- read.table("precip.txt", "\t", header=T)

# Import PSM data for 6 sites
psm <- read.table("psm.txt","\t",header=T)

# Merge predictor layer data with PSM data
psm <- merge(psm, catchments.final[,c("Model.catch","Impervious","dense.urban","light.med.urb","grass.shrub.crops","conifer.deciduous",
	"Apt.condo","Commercial","Industrial","Parks.open.space","Residential","FC00.Road","Nonlocal.Road")],
	by.x = "site", by.y = "Model.catch", all.x = T) 

psm <- merge(psm, precip, by = "year", all.x = T)
psm <- psm[order(psm$site, psm$year),]

# Import all catchments (for predictions)
catchments.predict <- read.table("catchments.predict.txt","\t",header=T)

# Import NMDS results (ordination axes) from PC-ORD output files
#UW2002.nms <- read.table("UW2002.nms.txt","\t",header=T)
#parcel.nms <- read.table("parcel.nms.txt","\t",header=T)
#roads.nms <- read.table("roads.nms.txt","\t",header=T)


#========================================================================================================================================
# ANALYSIS WITH FULL SUITE OF LAND USE / LAND COVER PREDICTORS
#========================================================================================================================================

#--------------------------------------------------------
# GENERALIZED LINEAR MIXED-EFFECTS MODELS FOR MORTALITY
#--------------------------------------------------------

# Response is binomial (n.psm/n.total)
# There is a random effect on the intercept, grouped by site
# Candidate models consist of combinations of at most 1 predictor from each of 4 categories
# (impervious, LU/LC, parcels, roads)

library(lme4)

# Create text strings of model terms
# Exclude any combinations of variables whose Spearman rank correlation is > 0.9 in absolute value
# conifer.deciduous, Impervious
# Residential, light.med.urb
# FC00.Road, light.med.urb
# FC00.Road, Parks.open.space
impervious.terms <- rep(c("", "Impervious"), each=90)
UW2002.terms <- rep(rep(c("","dense.urban","light.med.urb","grass.shrub.crops","conifer.deciduous"), each=18), 2)
parcel.terms <- rep(rep(c("","Apt.condo","Commercial","Industrial","Parks.open.space","Residential"), each=3), 10)
roads.terms <- rep(c("","FC00.Road","Nonlocal.Road"), 60)
 
model.terms <- paste(impervious.terms, UW2002.terms, parcel.terms, roads.terms, sep=" ")

for(i in 1:length(model.terms))
{
	blah <- unlist(strsplit(model.terms[i], " "))
	model.terms[i] <- ifelse((any(blah == "Impervious") & any(blah == "conifer.deciduous")) | (any(blah == "Residential") & any(blah == "light.med.urb")) |
				(any(blah == "FC00.Road") & any(blah == "light.med.urb")) | (any(blah == "FC00.Road") & any(blah == "Parks.open.space")),
				NA, paste(blah[blah != ""], collapse = " + "))
}
rm(blah)
model.terms[1] <- "1"
model.terms <- model.terms[!is.na(model.terms)]


# Use text strings to build list of models
psm.glmer.list <- vector("list",length(model.terms))

for(i in 1:length(psm.glmer.list))
	psm.glmer.list[[i]] <- glmer(formula(paste("psm ~ ", model.terms[i], " + (1 | site)")), data = psm, weights = n.total, family = binomial)

# Data frame of model comparisons
psm.glmer.tab <- data.frame(model=model.terms, k=NA, n=NA, Intercept=NA, RE.sd=NA, Impervious=NA, dense.urban=NA, light.med.urb=NA,
	grass.shrub.crops=NA, conifer.deciduous=NA, Apt.condo=NA, Commercial=NA, Industrial=NA, Parks.open.space=NA, Residential=NA,
	FC00.Road=NA, Nonlocal.Road=NA, deviance=NA, AICc=NA, dAICc=NA, wAICc=NA)
psm.glmer.tab$k <- sapply(psm.glmer.list, function(x) length(fixef(x)) + 1)
psm.glmer.tab$n <- sapply(psm.glmer.list, function(x) x@dims["n"])
psm.glmer.tab$Intercept <- sapply(psm.glmer.list, function(x) fixef(x)[1])
psm.glmer.tab$RE.sd <- sapply(psm.glmer.list, function(x) x@ST[[1]][1])
for(i in 6:17)
	psm.glmer.tab[,i] <- sapply(psm.glmer.list, function(x) fixef(x)[names(psm.glmer.tab)[i]])
psm.glmer.tab$deviance <- sapply(psm.glmer.list, deviance)
psm.glmer.tab$AICc <- psm.glmer.tab$deviance + 2*psm.glmer.tab$k + 2*psm.glmer.tab$k*(psm.glmer.tab$k + 1)/(psm.glmer.tab$n - psm.glmer.tab$k - 1)
psm.glmer.tab$dAICc <- psm.glmer.tab$AICc - min(psm.glmer.tab$AICc)
psm.glmer.tab$wAICc <- exp(-psm.glmer.tab$dAICc/2)/sum(exp(-psm.glmer.tab$dAICc/2))
psm.glmer.tab <- psm.glmer.tab[order(psm.glmer.tab$dAICc),]

# Calculate variable weights across 95% confidence set of models
mod95 <- psm.glmer.tab[cumsum(psm.glmer.tab$wAICc) < 0.95, ]
mod95$wAICc <- mod95$wAICc/sum(mod95$wAICc)
psm.var.wts <- t(as.matrix(mod95$wAICc))%*%ifelse(is.na(mod95[,6:17]), 0, 1)
psm.var.wts <- data.frame(wAICc=t(psm.var.wts))

# Calculate model-averaged fixed effect coefficients
psm.coef.ma <- as.matrix(mod95[,6:17])
psm.coef.ma[is.na(psm.coef.ma)] <- 0
psm.coef.ma <- t(as.matrix(mod95$wAICc))%*%psm.coef.ma

# Data frame with model avg estimates, SEs and var wts
psm.coefs <- data.frame(psm.var.wts, modavg.est = as.vector(psm.coef.ma), modavg.SE = NA)

# Unconditional SEs
var.theta <- lapply(psm.glmer.list[as.numeric(row.names(mod95))], function(x) summary(x)@coefs)
for(i in 1:nrow(psm.coefs))
{
 var.theta.i <- sapply(var.theta, function(x) ifelse(any(row.names(x) == row.names(psm.coefs)[i]), x[row.names(psm.coefs)[i], "Std. Error"], NA))
 bias.squared.i <- (sapply(psm.glmer.list[as.numeric(row.names(mod95))], function(x) fixef(x)[row.names(psm.coefs)[i]]) - psm.coefs$modavg.est[i])^2
 wts.i <- mod95$wAICc[!is.na(var.theta.i)]/sum(mod95$wAICc[!is.na(var.theta.i)])
 psm.coefs$modavg.SE[i] <- sum(sqrt(var.theta.i[!is.na(var.theta.i)] + bias.squared.i[!is.na(bias.squared.i)])*wts.i)
 rm(list=c("var.theta.i","bias.squared.i","wts.i"))
}
rm(var.theta)


# "Predictive map"
# Calculate model-averaged predicted PSM for unsampled catchments.
# These predictions are for the average population (intercept random effect = 0),
# so in a sense they are at the metapopulation or ESU level, not specific to any
# particular popn.

psm.predict <- data.frame(catchments.predict, psm.predicted.ma=0)
newdat <- as.matrix(catchments.predict[,match(row.names(psm.coefs), names(catchments.predict))])
for(i in 1:nrow(mod95))
{
 coefs.i <- psm.glmer.tab[i,match(row.names(psm.coefs), names(psm.glmer.tab))]
 coefs.i[is.na(coefs.i)] <- 0
 logit.psm <- as.vector(psm.glmer.tab$Intercept[i] + newdat%*%t(coefs.i))
 psm.i <- exp(logit.psm)/(1 + exp(logit.psm))
 psm.predict$psm.predicted.ma <- psm.predict$psm.predicted.ma + mod95$wAICc[i]*psm.i
 rm(list=c("coefs.i","logit.psm","psm.i"))
}
rm(newdat)
write.table(psm.predict, "psm.predict.txt", sep="\t")



#========================================================================================================================================
# ANALYSIS WITH RESTRICTED LAND USE / LAND COVER PREDICTORS AND PRECIPITATION
#========================================================================================================================================

#--------------------------------------------------------
# GENERALIZED LINEAR MIXED-EFFECTS MODELS FOR MORTALITY
#--------------------------------------------------------

# Response is binomial (n.psm/n.total)
# There is a random effect on the intercept, grouped by site
# Candidate models consist of combinations of the predictors:
# Impervious, FC00.Road, dry, rain, Impervious:dry, Impervious:rain, FC00.Road:dry, FC00.Road:rain
# with the restriction that Impervious OR FC00.Road may appear in a model, but not both.
# Premise: compare hypotheses that each of these factors is the best predictor of PSM.

library(lme4)

# Create text strings of model terms
model.terms2 <- c("1", "Impervious", "FC00.Road", "dry", "rain",
"Impervious + dry", "Impervious + rain", "FC00.Road + dry", "FC00.Road + rain", "dry + rain",
"Impervious + dry + Impervious:dry", "Impervious + rain + Impervious:rain", "FC00.Road + dry + FC00.Road:dry", "FC00.Road + rain + FC00.Road:rain",

"Impervious + dry + rain", "FC00.Road + dry + rain",
"Impervious + dry + rain + Impervious:dry", "Impervious + dry + rain + Impervious:rain", "Impervious + dry + rain + Impervious:dry + Impervious:rain",
"FC00.Road + dry + rain + FC00.Road:dry", "FC00.Road + dry + rain + FC00.Road:rain", "FC00.Road + dry + rain + FC00.Road:dry + FC00.Road:rain")


# Use text strings to build list of models
psm.glmer.list2 <- vector("list",length(model.terms2))

for(i in 1:length(psm.glmer.list2))
	psm.glmer.list2[[i]] <- glmer(formula(paste("psm ~ ", model.terms2[i], " + (1 | site)")), data = psm, weights = n.total, family = binomial)

# Data frame of model comparisons
psm.glmer.tab2 <- data.frame(model=model.terms2, k=NA, n=NA, Intercept=NA, RE.sd=NA, Impervious=NA, FC00.Road=NA, dry=NA, rain=NA,
	ImperviousXdry=NA, ImperviousXrain=NA, FC00.RoadXdry=NA, FC00.RoadXrain=NA, deviance=NA, AICc=NA, dAICc=NA, wAICc=NA)
names(psm.glmer.tab2)[10:13] <- c("Impervious:dry", "Impervious:rain", "FC00.Road:dry", "FC00.Road:rain")
psm.glmer.tab2$k <- sapply(psm.glmer.list2, function(x) length(fixef(x)) + 1)
psm.glmer.tab2$n <- sapply(psm.glmer.list2, function(x) x@dims["n"])
psm.glmer.tab2$Intercept <- sapply(psm.glmer.list2, function(x) fixef(x)[1])
psm.glmer.tab2$RE.sd <- sapply(psm.glmer.list2, function(x) x@ST[[1]][1])
for(i in 6:13)
	psm.glmer.tab2[,i] <- sapply(psm.glmer.list2, function(x) fixef(x)[names(psm.glmer.tab2)[i]])
psm.glmer.tab2$deviance <- sapply(psm.glmer.list2, deviance)
psm.glmer.tab2$AICc <- psm.glmer.tab2$deviance + 2*psm.glmer.tab2$k + 2*psm.glmer.tab2$k*(psm.glmer.tab2$k + 1)/(psm.glmer.tab2$n - psm.glmer.tab2$k - 1)
psm.glmer.tab2$dAICc <- psm.glmer.tab2$AICc - min(psm.glmer.tab2$AICc)
psm.glmer.tab2$wAICc <- exp(-psm.glmer.tab2$dAICc/2)/sum(exp(-psm.glmer.tab2$dAICc/2))
psm.glmer.tab2 <- psm.glmer.tab2[order(psm.glmer.tab2$dAICc),]


# Calculate variable weights across 95% confidence set of models
mod95.2 <- psm.glmer.tab2[cumsum(psm.glmer.tab2$wAICc) < 0.95, ]
mod95.2$wAICc <- mod95.2$wAICc/sum(mod95.2$wAICc)
psm.var.wts2 <- t(as.matrix(mod95.2$wAICc))%*%ifelse(is.na(mod95.2[,6:13]), 0, 1)
psm.var.wts2 <- data.frame(wAICc=t(psm.var.wts2))

# Calculate model-averaged fixed effect coefficients
psm.coef.ma2 <- as.matrix(mod95.2[,6:13])
psm.coef.ma2[is.na(psm.coef.ma2)] <- 0
psm.coef.ma2 <- t(as.matrix(mod95.2$wAICc))%*%psm.coef.ma2

# Data frame with model avg estimates, SEs and var wts
psm.coefs2 <- data.frame(psm.var.wts2, modavg.est = as.vector(psm.coef.ma2), modavg.SE = NA)

# Unconditional SEs
var.theta2 <- lapply(psm.glmer.list2[as.numeric(row.names(mod95.2))], function(x) summary(x)@coefs)
for(i in 1:nrow(psm.coefs2))
{
 var.theta.i <- sapply(var.theta2, function(x) ifelse(any(row.names(x) == row.names(psm.coefs2)[i]), x[row.names(psm.coefs2)[i], "Std. Error"]^2, NA))
 wts.i <- mod95.2$wAICc[!is.na(var.theta.i)]/sum(mod95.2$wAICc[!is.na(var.theta.i)])
 ma.est.i <- mod95.2[,row.names(psm.coefs2)[i]]*mod95.2$wAICc
 ma.est.i <- sum(ma.est.i, na.rm=T)/sum(mod95.2$wAICc[!is.na(ma.est.i)])	# model-avg estimate, conditional on the param being in the model
 bias.squared.i <- (mod95.2[,row.names(psm.coefs2)[i]] - ma.est.i)^2
 psm.coefs2$modavg.SE[i] <- sum(sqrt(var.theta.i[!is.na(var.theta.i)] + bias.squared.i[!is.na(bias.squared.i)])*wts.i)
 rm(list=c("var.theta.i","ma.est.i","bias.squared.i","wts.i"))
}
rm(var.theta2)


#--------------------------------------------------------
# GRAPHS
#--------------------------------------------------------

# Plot of observed and fitted PSM vs dry-season precip, using symbols to indicate sites
plot(psm$dry, psm$psm, xlab="Jun-Aug rainfall (in)", ylab="Prespawn mortality", pch="")
sym <- c(0,1,2,3,6,8)
cols <- c("black","blue","green","red","orange","purple")
for(i in 1:length(levels(psm$site)))
{
 points(psm$dry[psm$site==levels(psm$site)[i]], psm$psm[psm$site==levels(psm$site)[i]], pch=sym[i], col=cols[i], cex=1.5)
 dry.range <- seq(min(psm$dry), max(psm$dry), length=100)
 FC00.Road.i <- psm$FC00.Road[psm$site==levels(psm$site)[i]][1]
 coef.i <- unlist(coef(psm.glmer.list2[[13]])$site[levels(psm$site)[i], ])
 logitfit <- coef.i["(Intercept)"] + FC00.Road.i*coef.i["FC00.Road"] + dry.range*coef.i["dry"] + FC00.Road.i*dry.range*coef.i["FC00.Road:dry"]
 fit <- exp(logitfit)/(1 + exp(logitfit))
 lines(dry.range, fit, col = cols[i], lwd = 2)
}





#------------------------------------------------------
# NONMETRIC MULTIDIMENSIONAL SCALING
#------------------------------------------------------

# Levels of Model.catch that correspond to study sites
study.sites <- c("Des Moines","Fauntleroy","Fortson","Longfellow","Pipers","Thornton")

# Plot ordinations

# LU/LC
dev.new()
par(mfrow=c(1,1))
plot(UW2002.nms$axis1, UW2002.nms$axis2, xlab="NMS Axis 1", ylab="NMS Axis 2", main = "Land use/land cover", cex=1.5, cex.lab=1.5, cex.axis=1.2)
points(UW2002.nms$axis1[is.element(UW2002.nms$Model.catch, study.sites)], UW2002.nms$axis2[is.element(UW2002.nms$Model.catch, study.sites)],
	pch=16, cex=1.5, col="red")

# Parcels
dev.new()
par(mfrow=c(1,1))
plot(parcel.nms$axis1, parcel.nms$axis2, xlab="NMS Axis 1", ylab="NMS Axis 2", main = "Parcel classifications", cex=1.5, cex.lab=1.5, cex.axis=1.2)
points(parcel.nms$axis1[is.element(parcel.nms$Model.catch, study.sites)], parcel.nms$axis2[is.element(parcel.nms$Model.catch, study.sites)],
	pch=16, cex=1.5, col="red")

# Roads
dev.new()
par(mfrow=c(1,1))
plot(roads.nms$axis1, roads.nms$axis2, xlab="NMS Axis 1", ylab="NMS Axis 2", main = "Road density", cex=1.5, cex.lab=1.5, cex.axis=1.2)
points(roads.nms$axis1[is.element(roads.nms$Model.catch, study.sites)], roads.nms$axis2[is.element(roads.nms$Model.catch, study.sites)],
	pch=16, cex=1.5, col="red")


# Plot NMDS axes against raw variables and examine rank correlations

# LU/LC
dev.new()
par(mfcol=c(2,6))
for(i in names(catchments.final)[18:23])
for(j in c("axis1","axis2"))
{
rij <- cor.test(catchments.final[,i], UW2002.nms[,j], method="spearman")
plot(asin(sqrt(catchments.final[,i])), UW2002.nms[,j], cex=1.5, xlab=paste(c("asin(sqrt(",i,"))"), collapse=""), ylab=j, 
	main=paste(c("rho = ", round(rij$est,2), " p = ", round(rij$p.val,4)), collapse=""))
}

# Parcels
dev.new()
par(mfcol=c(2,6))
for(i in names(catchments.final)[24:29])
for(j in c("axis1","axis2"))
{
rij <- cor.test(catchments.final[,i], parcel.nms[,j], method="spearman")
plot(asin(sqrt(catchments.final[,i])), parcel.nms[,j], cex=1.5, xlab=paste(c("asin(sqrt(",i,"))"), collapse=""), ylab=j, 
	main=paste(c("rho = ", round(rij$est,2), " p = ", round(rij$p.val,4)), collapse=""))
}

# Roads
dev.new()
par(mfcol=c(2,7))
for(i in names(catchments.final)[30:36])
for(j in c("axis1","axis2"))
{
rij <- cor.test(catchments.final[,i], roads.nms[,j], method="spearman")
plot(asin(sqrt(catchments.final[,i])), roads.nms[,j], cex=1.5, xlab=paste(c("asin(sqrt(",i,"))"), collapse=""), ylab=j, 
	main=paste(c("rho = ", round(rij$est,2), " p = ", round(rij$p.val,4)), collapse=""))
}






