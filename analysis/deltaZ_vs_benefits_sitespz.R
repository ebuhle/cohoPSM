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
## ailene.ettinger@tnc.org
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
#set_here("/Users/aileneettinger/Documents/GitHub/cohoPSM")

## Step 1:  Read in the data/model estimates- use only predicted attributes for now
psm_pre <- read.table(here("analysis","results","PSM_predictions.txt"), header=TRUE)
spawn <- read.csv(here("data","spawner_data.csv"), header=TRUE)
spatial_pred <- read.csv(here("data","spatial_data_predict.csv"), header=TRUE)
z<-read.csv(here("analysis","results","psm_z_threshold_colors.csv"), header=TRUE)
#spatial<- read.csv(here("data","spatial_data.csv"), header=TRUE)

salmon <- read.csv(here("data","salmonscape","WA_integrated_Fish_coho_chinook_chum_with_WADOE.csv"),header=TRUE, skip=1)

## Step 2: Choices: select the threshold psm and you want to use, and select all IDs or only IDs for which we have PSM data (rather than predicted PSM)
input <- as.data.frame(NA)
input$psm_thresh <-psm_thresh<- 0.3
input$attribute<-"Coho_Presence_m"
predsites <- TRUE #if false, selects out only IDs with PSM calculated from field data, rather than IDs with predicted PSM

## Step 3: combine all the data and prep for plotting calculate mean spawner abundance by ID, across years
## combined data file with things we want to plot is called "d"
source(here("analysis","source","prepforplots.R"))

dim(d)

## Step 4. Function to plot change in Z on x-axis and attribute of interest on the y axis
zplotfx <- function(psm_thresh,attribut){
  figname<-paste(attribut,psm_thresh,".pdf",sep="_")
  #Commented out my way of calculating deltaz, as using Eric's for now
  #Zcrit<-min(d$Z[d$p_psm_mean>psm_thresh], na.rm=TRUE)
  #Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
  # for all sites and select out just the bad sites
  d$Zcrit<-Zcrit
  d$deltaZ<-d$delta_Z
   #add a column for the colors to plot for whether or not site has below threshold psm
  d$psmshape<-17 
  d$psmshape[d$p_psm_mean<psm_thresh]<-19
  attribute<-d[,which(colnames(d)==attribut)]
  d$attribut_stan<-(attribute-mean(attribute,na.rm=TRUE))/sd(attribute,na.rm=TRUE)
  dxy<-subset(d,select=c(delta_Z,attribut_stan))
  score<-as.matrix(dist(rbind(c(0,max(dxy$attribut_stan,na.rm=TRUE)),dxy), method="euclidean"))[1,-1]
  dxy<-cbind(d$ID,dxy,d[,which(colnames(d)==attribut)],score)
  dxy<-dxy[-which(is.na(dxy$attribut_stan)),]
  dxy<-dxy[order(dxy$score),]
  restPalette <- colorRampPalette(brewer.pal(9, "Blues")) #### Gives us a heat map look
  consPalette <- colorRampPalette(brewer.pal(9, "Greens")) #### Gives us a heat map look
  dxy$cols<-NA
  dxy$cols[which(dxy$delta_Z<0)]<-rev(restPalette(length(dxy$score[which(dxy$delta_Z<=0)])))
  dxy$cols[which(dxy$delta_Z>0)] <- rev(consPalette(length(dxy$score[which(dxy$delta_Z>0)])))
  
  #dxy<- data.frame(cbind(dxy,cols))
  colnames(dxy)[1:4]<-c("ID","delta_Z","benefit.stan","benefit")
  #quartz(width = 8, height = 5)
  pdf(here("analysis","results","figures",figname))
  plot(dxy$delta_Z,dxy$benefit, cex=1.5,cex.lab=1.2,cex.axis=1.2,xlab="Effort (Delta Z)", ylab= paste(attribut), type="p", pch=d$psmshape, col=dxy$cols)
  abline(v=0, lty=2, lwd=2, col="black")
  #text(Zcrit,max(dxy$benefit),"zcrit", col= "blue")
  #mtext(side=1,"high",line=4,adj=1,cex=0.8)
  #mtext(side=1,"low",line=4,adj=0,cex=0.8)
  write.csv(dxy,here("analysis","results",paste(attribut,"scores.csv", sep="_")), row.names = FALSE)
  #legend("topleft", legend=c("Highest priority","Lowest priority"), pch=19,col=c(cols[1],cols[length(cols)]), cex=.8, bty="n")
  dev.off()
  }


zplotfx (input$psm_thresh,"Coho_Presence_m")
zplotfx (input$psm_thresh,"nsp_pres")
zplotfx (input$psm_thresh,"ChinFa_Presence_m")


##Make a schematic diagram showing our approach 
#dev.new(height=8,width=16)
#pdf(here("analysis","results","figures","benefits_stream_spp_wscheme.pdf"), width = 15, height = 8)
pdf(here("analysis","results","figures","scheme.pdf"), width = 8, height = 5)
psm_thresh<-0.30
Zcrit<-min(d$Z_mean[d$p_psm_mean>psm_thresh], na.rm=TRUE)
#quartz(height=7,width=14)
#par(mfrow=c(2,2))
#make a blank plot with the prioritization scheme
plot(psm_pre$deln[1:51],d$p_psm_mean[1:51], pch=19, col="white", yaxt='n',xaxt='n',cex.lab=1.2,cex.axis=1.5,cex=1.52, 
     xlab="Urbanization (Z)", ylab= "Biological attribute of interest")
mtext(side=1,"high",line=3,adj=1,cex=0.8)
mtext(side=1,"low",line=3,adj=0,cex=0.8)
mtext(side=3,"Restoration Priority",line=-4,adj=.9,cex=.9)
mtext(side=3,"Conservation Priority",line=-4,adj=.1,cex=.9)
mtext(side=1,"No action needed",line=-4,adj=.1,cex=.9)
mtext(side=1,"Low Ecological Priority",line=-4,adj=.9,cex=.9)
abline(v=0.3, lty=2, lwd=2, col="blue")
text(Zcrit+.04,.02,label="Zcrit", col="blue",cex=1.2)

polygon(psm_pre$Z_mean[1:51])
dev.off()
#plot relationship of PSM and Z
pdf(here("analysis","results","figures","psmvsZ.pdf"), width = 8, height = 5)
quartz()
plot(psm_pre$Z_mean[1:51],psm_pre$p_psm_mean[1:51],pch=19, col="gray", cex.lab=1.2,cex.axis=1.2,cex=1.52, 
     xlab="Urbanization score (Z)", ylab= "Pre-Spawn Mortality")

abline(h=input$psm_thresh, lty=2, lwd=2)
text(min(d$Z_mean, na.rm=TRUE)+1,input$psm_thresh+.02,label="PSM threshold", cex=1.2)
abline(v=Zcrit, lty=2, lwd=2, col="blue")
#text(psm_pre3$p.psm.mean, psm_pre3$Z.mean, labels=as.numeric(as.factor(psm_pre3$ID)),cex=0.8, font=2)
# polygon(c(Zcrit,Zcrit,max(d$Z_mean, na.rm=TRUE)+.5,max(d$Z_mean, na.rm=TRUE)+.5),c(input$psm_thresh,1,1,input$psm_thresh),
#         col=adjustcolor("salmon",alpha.f=0.5),
#         border=NA)
text(Zcrit+.04,.02,label="Zcrit", col="blue",cex=1.2)
mtext(side=1,"high",line=3,adj=1,cex=0.8)
mtext(side=1,"low",line=3,adj=0,cex=0.8)
#mtext(side=3,"Restoration",line=0,adj=1,cex=0.8)
#mtext(side=3,"Conservation",line=0,adj=0,cex=0.8)

dev.off()

#write a csv files with all benefits and their scores and color assignments to send to blake
Coho_Presence_m_scores
score_chinsppres_m<-read.csv(here("analysis","results","ChinFa_Presence_m_scores.csv"), header=TRUE)
score_cohosppres_m<-read.csv(here("analysis","results","Coho_Presence_m_scores.csv"), header=TRUE)
score_nsp<-read.csv(here("analysis","results","nsp_pres_scores.csv"), header=TRUE)

head(score_chinsppres_m)
head(score_nsp)
head(score_cohosppres_m)
score_chinsppres_m<-dxy
colnames(score_chinsppres_m)[3:6]<-c("chinsppres.stan","chinsp.pres.m","score.chinsp.pres","cols.chinsppres")
colnames(score_cohosppres_m)[3:6]<-c("coho.pres.m","coho.pres.stan","score.coho.pres","cols.cohosppres")

scores<-full_join(score_cohosppres_m,score_chinsppres_m)
colnames(scores)[3:5]<-c("coho.pres.m","coho.pres.stan","score.coho.pres")


#fix it up for Blake: remove "NAs" and replace with -9999 and remove colums that we don't care about
scores_forblake<-subset(scores,select=c(ID,Z,coho.pres.m,score.coho.pres,chinfa.pres.m,score.chinfa.pres,chinsp.pres.m,score.chinsp.pres,numspp,score.numsp))
scores_forblake[is.na(scores_forblake)]<-"-9999"

write.csv(scores_forblake,file=here("analysis","results","scores.csv"), row.names = FALSE)

