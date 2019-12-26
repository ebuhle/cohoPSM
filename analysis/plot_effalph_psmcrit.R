#make a table andplot showing how change sin alpha and psmcrit affect prioritization
#started by ailene
#dec 22, 2019

## housekeeping
rm(list=ls()) 
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
options(stringsAsFactors = FALSE)
##set working directory
setwd("~/GitHub/cohoPSM")
#setwd("/Users/aileneettinger/Documents/GitHub/cohoPSM")
## libraries
library(here)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
library(scales)
library(tidyr)
## Step 1: Make a table of the top highest prioritiy site for conservation and restoration, using different priorities

#Read in all the files in the "scores" folder and put together into one dataframe
files<-list.files(here("analysis","results","scores"))
#select out only the coho scores first
files<-files[grep("score.csv",files)]

allfiles<-c()
top25res<-c()
top25cons<-c()
alphpsmcrit<-c()
for(i in 1:length(files)){
f<-read.csv(here("analysis","results","scores",files[i]))
f<-f[order(f$score),]
f$alpha<-substr(files[i],1,4)
f$psmcrit<-substr(files[i],nchar(files[i])-12,nchar(files[i])-10)
allfiles<-rbind(allfiles,f)
frest<-f[f$delta_Z<0,]
frest<-frest[order(frest$score),]
frest$rank<-order(frest$score)
fcons<-f[f$delta_Z>0,]
fcons<-fcons[order(fcons$score),]
fcons$rank<-order(fcons$score)
top25res<-cbind(top25res,frest$ID[1:25])
top25cons<-cbind(top25cons,fcons$ID[1:25])

alphpsmcrit<-c(alphpsmcrit,substr(files[i],1,8))
}
colnames(top25res)<-colnames(top25cons)<-alphpsmcrit
rownames(top25res)<-rownames(top25cons)<-NULL

#convert allfiles from long format to wide format
allfiles<-allfiles[,-which(colnames(allfiles)=="cols")]
allfiles$ID<-factor(allfiles$ID)
allfiles$alpha<-factor(allfiles$alpha)

allfiles$psmcrit<-factor(allfiles$psmcrit)
allfiles<-allfiles[order(allfiles$ID,allfiles$alpha,allfiles$psmcrit),]
allfiles$action<-"cons"
allfiles$action[allfiles$delta_Z<0]<-"rest"
allsub<-subset(allfiles,select=c(ID,score,alpha,psmcrit))
allwide<-spread(allsub,alpha,score)
wide.3<-allwide[allwide$psmcrit==0.3,]
wide.4<-allwide[allwide$psmcrit==0.4,]
getdiffs<-function(allwide,psmcrit){
  w<-allwide[allwide$psmcrit==psmcrit,]
  w<-w[order(w$ID),]
  w$diff.9<-w$`0.95`-w$`0.9_`
  w$diff.8<-w$`0.95`-w$`0.8_`
 return(w) 
}
wide.2<-getdiffs(allwide,0.2)
wide.3<-getdiffs(allwide,0.3)
wide.4<-getdiffs(allwide,0.4)
wide.3$diff.alpha.4<-wide.3$`0.95`-wide.4$`0.95`
wide.3$diff.alpha.2<-wide.3$`0.95`-wide.2$`0.95`

#read in the score using chin as metric
chin<-read.csv(here("analysis","results","scores","0.95_0.3_scorechin.csv"))
chinsub<-subset(chin,select=c(ID,score))
chinsub$ID<-as.integer(chinsub$ID)
nsp<-read.csv(here("analysis","results","scores","0.95_0.3_scorensp.csv"))
nspsub<-subset(nsp,select=c(ID,score))
wide.3$ID<-as.integer(wide.3$ID)
chincoho<-left_join(wide.3,chinsub)
chincoho<-chincoho[-which(is.na(chincoho$score)),]
nspcoho<-left_join(wide.3,nspsub)
nspcoho<-nspcoho[-which(is.na(nspcoho$score)),]

windows(width=10,height=15)
par(mfrow=c(3,2))
plot(wide.3$`0.95`,wide.3$`0.9_`)
plot(wide.3$`0.95`,wide.3$`0.8_`)
plot(wide.3$`0.95`,wide.2$`0.95`)
plot(wide.3$`0.95`,wide.4$`0.95`)
plot(chincoho$`0.95`,chincoho$score)
plot(nspcoho$`0.95`,nspcoho$score)

#plot diffs on a single plot

