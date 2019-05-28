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
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(stringsAsFactors = FALSE)

## libraries
library(here)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
## Step 1:  Read in the data/model estimates
psm_pre <- read.table(here("analysis","results","PSM_predictions.txt"), header=TRUE)
spawn <- read.csv(here("data","spawner_data.csv"), header=TRUE)
spatial_pred <- read.csv(here("data","spatial_data_predict.csv"), header=TRUE)
spatial<- read.csv(here("data","spatial_data.csv"), header=TRUE)

salmon <- read.csv(here("data","salmonscape","WA_integrated_Fish_dist_with_WADOE_basins.csv"),header=TRUE)

## Step 2: Choices: select the threshold psm and you want to use, and select all sites or only sites for which we have PSM data (rather than predicted PSM)
input <- as.data.frame(NA)
input$psm_thresh <- 0.25
allsites <- TRUE #if false, selects out only sites with PSM calculated from field data, rather than sites with predicted PSM too

## Step 3: combine all the data and prep for plotting calculate mean spawner abundance by site, across years
## combined data file with things we want to plot is called "d"
source(here("analysis","source","prepforplots.R"))


#source ("../analysis/source/prepforplots.R")

dim(d)
## Step 4. Plot Change in Z on x-axis and benefits of interest on the y axis

#dev.new(height=8,width=16)
pdf(here("analysis","results","testdeltaZvsbenefitsfig.pdf"), width = 16, height = 8)
quartz()
par(mfrow=c(1,3))
#plot relationship of PSM and Z

plot(d$p_psm_mean, d$Z_mean, pch=19,col=d$psmcol, cex.lab=1.2,cex.axis=1.2,cex=2, 
     ylab="Urbanization score (Z)", xlab= "Mean Pred.PSM")
abline(v=input$psm_thresh, lty=2, lwd=2)
text(input$psm_thresh+.02,min(d$Z_mean),label="PSM threshold", cex=1.2)
abline(h=Zcrit, lty=2, lwd=2, col="blue")
#text(psm_pre3$p.psm.mean, psm_pre3$Z.mean, labels=as.numeric(as.factor(psm_pre3$site)),cex=0.8, font=2)
polygon(c(input$psm_thresh,1,1,input$psm_thresh),c(Zcrit,Zcrit,max(d$Z_mean)+.5,max(d$Z_mean)+.5),
        col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
text(.02,Zcrit+.04,label="Zcrit", col="blue",cex=1.2)

#plot change in Z vs benefit (first  benefit=spawner abundance)
plot(d$Z_mean, d$spawn.n, cex=2,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = spawner abundance", type="p", pch=19, col=d$psmcol)
     ## Step 4. Compare to other benefits of interest- abundance or presence of other salmon species for all sites, human pops, stream biodiversity. Can blake get these?
text(d$Z_mean, d$spawn.n, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Z_mean,na.rm=TRUE),max(d$Z_mean,na.rm=TRUE)),
        c(Zcrit,max(d$spawn.n, na.rm=TRUE),max(d$spawn.n, na.rm=TRUE),Zcrit),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")


head(d)

#plot change in Z vs benefit of number of people living nearby (don't have this data right now, so making it up as a function of Z with some noist)
#plot(psm_pre3$Z.mean,pop.pred)#check
plot(d$Z_mean,d$pop_census, cex=2,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = number of people nearby", type="p", pch=19, col=d$psmcol)
text(d$Z_mean, d$pop_census, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Z_mean,na.rm=TRUE),max(d$Z_mean,na.rm=TRUE)),
        c(Zcrit,max(d$pop_census, na.rm=TRUE),max(d$pop_census, na.rm=TRUE),Zcrit),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")

mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)

quartz()
par(mfrow=c(1,3))


#plot meters of stream with coho presence
plot(d$Z_mean,d$pop_census, cex=2,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = number of people nearby", type="p", pch=19, col=d$psmcol)
text(d$Z_mean, d$pop_census, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Z_mean),max(d$Z_mean)),
        c(Zcrit,max(d$pop_census, na.rm=TRUE),max(d$pop_census, na.rm=TRUE),Zcrit),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)

mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)
abline(v=Zcrit, lty=2, lwd=2, col="blue")




quartz()
par(mfrow=c(1,3))


#plot meters of stream with coho present
plot(d$Z_mean,d$Presence..m., cex=2,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = m stream with coho present", type="p", pch=19, col=d$psmcol)
#text(d$Z_mean, d$Presence..m, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Z_mean, na.rm=TRUE),max(d$Z_mean, na.rm=TRUE)),
        c(Zcrit,max(d$Presence..m, na.rm=TRUE),max(d$Presence..m, na.rm=TRUE),Zcrit),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")

mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)


#plot meters of stream with rearing coho
plot(d$Z_mean,d$Coho_Rear,cex=2,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = m stream with coho rearing", type="p", pch=19, col=d$psmcol)
#text(d$Z_mean, d$Rearing..m, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Z_mean, na.rm=TRUE),max(d$Z_mean, na.rm=TRUE)),
        c(Zcrit,max(d$Coho_Rear, na.rm=TRUE),max(d$Coho_Rear, na.rm=TRUE),Zcrit),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")

mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)

#plot meters of stream with spawning coho
plot(d$Z_mean,d$Spawning..m.,cex=2,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = m stream with coho spawning", type="p", pch=19, col=d$psmcol)
#text(d$Z_mean, d$Spawning..m, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Spawning..m., na.rm=TRUE),max(d$Spawning..m., na.rm=TRUE)),
        c(Zcrit,max(d$Spawning..m., na.rm=TRUE),max(d$Spawning..m., na.rm=TRUE),Zcrit),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")

mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)






dev.off()



