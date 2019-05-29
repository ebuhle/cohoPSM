## Estimate reduction in urbanization ("Z") required to reach 40% PSM for a particular ID.
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
## Step 1:  Read in the data/model estimates- use only predicted attributes for now
psm_pre <- read.table(here("analysis","results","PSM_predictions.txt"), header=TRUE)
spawn <- read.csv(here("data","spawner_data.csv"), header=TRUE)
spatial_pred <- read.csv(here("data","spatial_data_predict.csv"), header=TRUE)
#spatial<- read.csv(here("data","spatial_data.csv"), header=TRUE)

salmon <- read.csv(here("data","salmonscape","WA_integrated_Fish_dist_with_WADOE_basins.csv"),header=TRUE)

## Step 2: Choices: select the threshold psm and you want to use, and select all IDs or only IDs for which we have PSM data (rather than predicted PSM)
input <- as.data.frame(NA)
input$psm_thresh <- 0.25
predsites <- TRUE #if false, selects out only IDs with PSM calculated from field data, rather than IDs with predicted PSM too

## Step 3: combine all the data and prep for plotting calculate mean spawner abundance by ID, across years
## combined data file with things we want to plot is called "d"
source(here("analysis","source","prepforplots.R"))


#source ("../analysis/source/prepforplots.R")

dim(d)
## Step 4. Plot Change in Z on x-axis and benefits of interest on the y axis

#dev.new(height=8,width=16)
pdf(here("analysis","results","testdeltaZvsbenefitsfig.pdf"), width = 16, height = 8)
quartz(height=16,width=16)
par(mfrow=c(2,2))
#make a blank plot with the prioritization scheme
plot(psm_pre$Z_mean[1:51],d$p_psm_mean[1:51], pch=19, col="white", yaxt='n',xaxt='n',cex.lab=1.5,cex.axis=1.5,cex=1.52, 
     xlab="Urbanization", ylab= "Benefit/Biological attribute of interest")
mtext(side=1,"high",line=3,adj=1,cex=0.8)
mtext(side=1,"low",line=3,adj=0,cex=0.8)
mtext(side=3,"Restoration",line=-4,adj=.9,cex=1.2)
mtext(side=3,"Conservation",line=-4,adj=.1,cex=1.2)
mtext(side=1,"No action needed",line=-4,adj=.1,cex=1.2)
mtext(side=1,"Low Ecological Priority",line=-4,adj=.9,cex=1.2)



#plot relationship of PSM and Z

plot(psm_pre$Z_mean[1:51],psm_pre$p_psm_mean[1:51],pch=19, col="gray", cex.lab=1.2,cex.axis=1.2,cex=1.52, 
     xlab="Urbanization score (Z)", ylab= "Pre-Spawn Mortality")

abline(h=input$psm_thresh, lty=2, lwd=2)
text(min(d$Z_mean, na.rm=TRUE)+1,input$psm_thresh+.02,label="PSM threshold", cex=1.2)
abline(v=Zcrit, lty=2, lwd=2, col="blue")
#text(psm_pre3$p.psm.mean, psm_pre3$Z.mean, labels=as.numeric(as.factor(psm_pre3$ID)),cex=0.8, font=2)
polygon(c(Zcrit,Zcrit,max(d$Z_mean, na.rm=TRUE)+.5,max(d$Z_mean, na.rm=TRUE)+.5),c(input$psm_thresh,1,1,input$psm_thresh),
        col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
text(Zcrit+.04,.02,label="Zcrit", col="blue",cex=1.2)
mtext(side=1,"high",line=3,adj=1,cex=0.8)
mtext(side=1,"low",line=3,adj=0,cex=0.8)
#mtext(side=3,"Restoration",line=0,adj=1,cex=0.8)
#mtext(side=3,"Conservation",line=0,adj=0,cex=0.8)



#Plot meters of coho spawning and number of species present, number of species spawning
#quartz()
#par(mfrow=c(1,2))
#plot meters of stream with coho present
#Give each ID a score
#standardizing the effort and the benefit, so that they are equally weighted...we can decide if we want to weight things differently.

d$Coho_Pres_stan<-(d$Coho_Presence_m-mean(d$Coho_Presence_m, na.rm=TRUE))/sd(d$Coho_Presence_m, na.rm=TRUE)
dxy<-subset(d,select=c(Z_mean,Coho_Pres_stan))

score<-as.matrix(dist(rbind(c(0,max(dxy$Coho_Pres_stan,na.rm=TRUE)),dxy), method="euclidean"))[1,-1]

dxy<-cbind(d$ID,dxy,d$Coho_Presence_m,score)
dxy<-dxy[order(dxy$score),]
myPalette <- colorRampPalette(brewer.pal(9, "YlGnBu")) #### Gives us a heat map look
cols = myPalette(length(score))
dxy<- data.frame(cbind(dxy,cols))
colnames(dxy)[1:4]<-c("ID","Z","benefit.stan","benefit")


plot(dxy$Z,dxy$benefit, cex=1.5,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = m stream with coho present", type="p", pch=d$psmshape, col=dxy$cols)
#text(d$Z_mean, d$Presence..m, labels=as.numeric(as.factor(d$ID)),cex=0.8, font=2)
#polygon(c(Zcrit,Zcrit,max(dxy$Z, na.rm=TRUE),max(dxy$Z, na.rm=TRUE)),
#        c(min(dxy$benefit, na.rm=TRUE),max(dxy$benefit, na.rm=TRUE),max(dxy$benefit, na.rm=TRUE),min(dxy$benefit, na.rm=TRUE)),col=adjustcolor("salmon",alpha.f=0.5),
#        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")
#mtext(side=3,"Restoration",line=-2,adj=1,cex=0.8)
#mtext(side=3,"Conservation",line=-2,adj=0,cex=0.8)
mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)
score_cohopres_m<-dxy
legend("topleft", legend=c("Lowest priority","Highest priority"), pch=19,col=c(cols[length(cols)],cols[1]), cex=.8)



dxy<-subset(d,select=c(Z_mean,nsp_pres))

score<-as.matrix(dist(rbind(c(0,max(dxy$nsp_pres,na.rm=TRUE)),dxy), method="euclidean"))[1,-1]

dxy<-cbind(d$ID,dxy,d$nsp_pres,score)
dxy<-dxy[order(dxy$score),]
myPalette <- colorRampPalette(brewer.pal(9, "YlGnBu")) #### Gives us a heat map look
cols = myPalette(length(score))
dxy<- data.frame(cbind(dxy,cols))
colnames(dxy)[1:4]<-c("ID","Z","benefit.stan","benefit")


plot(dxy$Z,dxy$benefit, cex=1.5,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= "Benefit = # salmon species present", type="p", pch=d$psmshape, col=dxy$cols)
#text(d$Z_mean, d$Presence..m, labels=as.numeric(as.factor(d$ID)),cex=0.8, font=2)
#polygon(c(Zcrit,Zcrit,max(dxy$Z, na.rm=TRUE),max(dxy$Z, na.rm=TRUE)),
#        c(min(dxy$benefit, na.rm=TRUE),max(dxy$benefit, na.rm=TRUE),max(dxy$benefit, na.rm=TRUE),min(dxy$benefit, na.rm=TRUE)),col=adjustcolor("salmon",alpha.f=0.5),
#        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")
#mtext(side=3,"Restoration",line=1,adj=1,cex=0.8)
#mtext(side=3,"Conservation",line=0,adj=0,cex=0.8)
mtext(side=1,"high",line=4,adj=1,cex=0.8)
mtext(side=1,"low",line=4,adj=0,cex=0.8)
score_salmonsp<-dxy
colnames(score_salmonsp)[3:6]<-c("numspp.stan","numspp","score.numsp","cols.numsp")
scores<-full_join(score_cohopres_m,score_salmonsp)
scores<-scores[,-6]
scores<-scores[,-9]
colnames(scores)[3:5]<-c("coho.pres.m","coho.pres.stan","score.coho.pres")
#fix it up for Blake: remove "NAs" and replace with -9999 and remove colums that we don't care about
scores_forblake<-subset(scores,select=c(ID,Z,coho.pres.m,score.coho.pres,numspp,score.numsp))
scores_forblake[is.na(scores_forblake)]<-"-9999"

write.csv(scores,file="output/scores.csv", row.names = FALSE)
head(scores)
