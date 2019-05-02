## Estimate reduction in urbanization ("Z") required to reach 40% PSM for a particular site.
## Goal is to map change in Z required to reach threshold PSM and relate it to benefits of interest
## Potential benefits include:
## 1. Number of coho
## 2. Presence of spawning/rearing coho and/or other slamon species
## 2. Numer of other salmon species
## 3. Number of stream species (biodiversity)
## 4. Number of ESA-protected stocks
## 5. Number of people living 1ithin 1 km of the stream
## Started by Ailene 1 May 2019
## ailene.ettinger@noaa.gov
#################################################################
#################################################################
## housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

## libraries
library(dplyr)
## Set working directory. Add in your own path in an if statement for your file structure
if(length(grep("ailene", getwd())>0)) { setwd("~/Documents/GitHub/cohoPSM")}

## Step 1: Read in the data and model estimates
psm_pre<-read.table("analysis/results/PSM_predictions.txt",header=TRUE)
head(psm_pre)

plot(psm_pre$Z.mean,psm_pre$p.psm.mean,)
#A cheap way to get Zcrit associated with threshold PSM. 
#Need to make this more robust, use model/include full posterior distribution to get error, etc
psm_thresh<-0.25

Zcrit<-mean(psm_pre$Z.mean[round(psm_pre$p.psm.mean,digits=2)==psm_thresh])

#try using a lower Zcrit- the lower Z value associated with a psm above the threshold
#Zcrit<-min(psm_pre$Z.mean[psm_pre$p.psm.mean>psm_thresh])

## Step 2. Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
## for all sites and select out just the bad sites
psm_pre$Zcrit<-Zcrit
psm_pre$deltaZ<-psm_pre$Zcrit-psm_pre$Z.mean

## Step 3. Plot Change in Z on x-axis and benefits of interest on the y axis
## as a mock up of what we may want to do, use mean abundance of spawning salmon as the measure of abundance.
spawn<-read.csv("data/spawner_data.csv", header=TRUE)

## calculate mean spawner abundance by site, across years
spawnmn<-aggregate(spawn$n,list(spawn$site),mean)
colnames(spawnmn)<-c("site", "spawn.n")

#used pared down psm_pre file, with just named sites (what are numbered sites anyway?)
psm_pre2<-psm_pre[1:51,]
psm_pre3<-full_join(psm_pre2,spawnmn)
#create a variable of 0/1 for if PSM is less than 0.4
psm_pre3$psmcol<-"darkred" 
psm_pre3$psmcol[psm_pre3$p.psm.mean<psm_thresh]<-"lightgreen"

#quartz(height=8,width=16)
pdf("analysis/results/testdeltaZvsbenefitsfig.pdf", width = 16, height = 8)

par(mfrow=c(1,3))
#plot relationship of PSM and Z

plot(psm_pre3$Z.mean, psm_pre3$p.psm.mean, pch=19,col=psm_pre3$psmcol, cex.lab=1.2,cex.axis=1.2,cex=2, xlab="Urbanization score (Z)", ylab= "Mean Pred.PSM")
abline(h=psm_thresh, lty=2, lwd=2)
text(min(psm_pre3$Z.mean)+.5,psm_thresh+.02,label="PSM threshold", cex=1.2)
abline(v=Zcrit, lty=2, lwd=2, col="blue")
text(Zcrit,.02,label="Zcrit", col="blue",cex=1.2)
text(psm_pre3$Z.mean, psm_pre3$p.psm.mean, labels=as.numeric(as.factor(psm_pre3$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(psm_pre3$Z.mean),max(psm_pre3$Z.mean)),
        c(0,1,1,0),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)

#plot change in Z vs benefit (first  benefit=spawner abundance)
plot(psm_pre3$deltaZ, psm_pre3$spawn.n, cex=2,cex.lab=1.2,cex.axis=1.2,xlab=expression(paste("Change in urbanization required (",Delta,"Z)", sep="")), ylab= "Benefit = spawner abundance", type="p", pch=19, col=psm_pre3$psmcol)
     ## Step 4. Compare to other benefits of interest- abundance or presence of other salmon species for all sites, human pops, stream biodiversity. Can blake get these?
text(psm_pre3$deltaZ, psm_pre3$spawn.n, labels=as.numeric(as.factor(psm_pre3$site)),cex=0.8, font=2)
polygon(c(0,0,min(psm_pre3$deltaZ),min(psm_pre3$deltaZ)),
        c(0,max(psm_pre3$spawn.n),max(psm_pre3$spawn.n),0),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
## Step 4. Compare to other benefits of interest- abundance or presence of other salmon species for all sites, human pops, stream biodiversity. Can blake get these?

#plot change in Z vs benefit of number of people living nearby (don't have this data right now, so making it up as a function of Z with some noist)
pop.pred<-as.integer(3000+psm_pre3$Z.mean*1000) + as.integer(rnorm(length(psm_pre3$Z.mean),0,500))
#plot(psm_pre3$Z.mean,pop.pred)#check
plot(psm_pre3$deltaZ,pop.pred, cex=2,cex.lab=1.2,cex.axis=1.2,xlab=expression(paste("Change in urbanization required (",Delta,"Z)", sep="")), ylab= "Benefit = number of people nearby", type="p", pch=19, col=psm_pre3$psmcol)
text(psm_pre3$deltaZ, pop.pred, labels=as.numeric(as.factor(psm_pre3$site)),cex=0.8, font=2)
polygon(c(0,0,min(psm_pre3$deltaZ),min(psm_pre3$deltaZ)),
        c(0,max(pop.pred),max(pop.pred),0),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)

dev.off()
## I looked into salmon abunndance data from DWF stock inventory pops- no "bad sites" included. some good sites included...look into this more
badsites<-psm_pre[psm_pre$Z_mean>psm_thresh,]
sitenames<-unique(psm_pre$site)[1:51]
badsites<-badsites[1:24,]#just the sites with names

abund<-read.csv("data/WDFW-Salmonid_Stock_Inventory_Populations.csv", header=TRUE)
sort(unique(abund$Population.Name))
dim(badsites)
tail(badsites)
badsites$site
#Questions:
#What are the numbered sites?