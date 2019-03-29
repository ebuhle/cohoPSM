#======================================================
# COHO PRESPAWN MORTALITY ANALYSIS
# (Feist et al. 2010, Science, in press)
#
# Ordination of catchment LU/LC data, and mixed
# models of mortality using ordination scores as
# predictors.
#======================================================

#------------------------------------------------------
# DATA IMPORT AND MANIPULATION
#------------------------------------------------------

# Import main catchment dataset (with 6 study sites included)
#PCAData <- read.table("/Users/blakefeist/Ranalyses/PSM/PCAData.csv",",",header=T)
PCAData <- read.table("PCAData.txt","\t",header=T)

# Import data for catchments that overlap with 6 study sites
omit.catch <- read.table("mds_catch_leave_out.txt","\t",header=T)

# Toss out catchments that intersect 6 study sites
PCAData2 <- PCAData[!is.element(PCAData$Model.catch, omit.catch$Model.catch),]

#toss out small catchments
PCAData2 <- PCAData2[PCAData2$Catch.Area.ft2 > 1e6,]

# Write out reduced dataset to a text file
write.table(PCAData2, "PCAData2.txt", row.names=F, sep="\t")

# After doing this, we can load the reduced version directly
PCAData2 <- read.table("PCAData2.txt","\t",header=T)


#------------------------------------------------------
# EXPLORATORY PLOTS
#------------------------------------------------------

#cdf of catchment areas (ft2)+
plot.ecdf(PCAData$Catch.Area.ft2, xlim = c(0, 8e8))
abline(v = 6e6)

#cumulative sum of catchment areas (ft2)
plot(1:nrow(PCAData), cumsum(sort(PCAData$Catch.Area.ft2))/sum(PCAData$Catch.Area.ft2),type = "l")

#histogram of catchment areas (ft2)
hist(PCAData$Catch.Area.ft2, 50)
abline(v = 120000000, col = "red")

#scatterplots for UW2002
pairs(PCAData2[,18:23])

#arcsine transform for UW2002
pairs(asin(sqrt(PCAData2[, 18:23])))

#scatterplots for parcels
pairs(PCAData2[,23:29])

#arcsine transform for parcels
pairs(asin(sqrt(PCAData2[,24:30])))

#scatterplots for roads
pairs(PCAData2[,31:40])

#arcsine transform for roads
pairs(asin(sqrt(PCAData2[,31:40])))

#------------------------------------------------------
# PCA
# (Don't use)
#------------------------------------------------------

# LU/LC
PCAUW2002 <- princomp(asin(sqrt(PCAData2[,18:23])),cor=T,scores=T)

# Parcels
PCAParcel <- princomp(asin(sqrt(PCAData2[,24:30])),cor=T,scores=T)

# Roads
PCARoad <- princomp(asin(sqrt(PCAData2[,31:40])),cor=T,scores=T)


#------------------------------------------------------
# MDS
#------------------------------------------------------

# Calculate Bray-Curtis dissimilarities for each suite of variables
#install.packages("ecodist")
#library(ecodist)

## This function in the ecodist package calculates B-C distances, but maxes out memory.
## Do it manually instead.
## distance(PCAData2[,18:23], "bray-curtis")

####################
# Create matrices
# ONLY NEED TO DO THIS ONCE!!!
bcmatUW2002 <- matrix(0, nrow = nrow(PCAData2), ncol = nrow(PCAData2))
bcmatParcel <- matrix(0, nrow = nrow(PCAData2), ncol = nrow(PCAData2))
bcmatRoad <- matrix(0, nrow = nrow(PCAData2), ncol = nrow(PCAData2))

# Loop over matrix elements and calculate pairwise distances
for(i in 1:nrow(PCAData2))
for(j in 1:i)
{
bcmatUW2002[i,j] <- sum(abs(PCAData2[i,18:23] - PCAData2[j,18:23]))/sum(PCAData2[i,18:23] + PCAData2[j,18:23])
bcmatUW2002[j,i] <- bcmatUW2002[i,j]

bcmatParcel[i,j] <- sum(abs(PCAData2[i,24:30] - PCAData2[j,24:30]))/sum(PCAData2[i,24:30] + PCAData2[j,24:30])
bcmatParcel[j,i] <- bcmatParcel[i,j]

bcmatRoad[i,j] <- sum(abs(PCAData2[i,31:40] - PCAData2[j,31:40]))/sum(PCAData2[i,31:40] + PCAData2[j,31:40])
bcmatRoad[j,i] <- bcmatRoad[i,j]
}
bcmatUW2002[is.nan(bcmatUW2002)] <- 0		# NaN occurs when two sites both have zeros for all variables. Dissimilarity should be 0.
bcmatParcel[is.nan(bcmatParcel)] <- 0		
bcmatRoad[is.nan(bcmatRoad)] <- 0

# Write out distance matrices
write(bcmatUW2002, "bcmatUW2002.txt", ncolumns=ncol(bcmatUW2002), sep="\t")
write(bcmatParcel, "bcmatParcel.txt", ncolumns=ncol(bcmatParcel), sep="\t")
write(bcmatRoad, "bcmatRoad.txt", ncolumns=ncol(bcmatRoad), sep="\t")
####################

# Convert matrices to dist objects (smaller representation of lower triangle)
bcdistUW2002 <- as.dist(bcmatUW2002)
bcdistParcel <- as.dist(bcmatParcel)
bcdistRoad <- as.dist(bcmatRoad)


# Do the ordinations
# Zero distances between distinct catchments are changed to a small number (1e-6)

# LU/LC
mdsUW2002 <- isoMDS(pmax(bcdistUW2002, 1e-6), k = 3, maxit = 500)

# Parcels
mdsParcel <- isoMDS(bcdistParcel, k = 3, maxit = 500)

# LU/LC
mdsRoad <- isoMDS(bcdistRoad, k = 3, maxit = 500)



