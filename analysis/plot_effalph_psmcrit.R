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
  w$diff.9<-abs((w$`0.95`-w$`0.9_`)/w$`0.95`)
  w$diff.8<-abs((w$`0.95`-w$`0.8_`)/w$`0.95`)
  return(w) 
}
wide.2<-getdiffs(allwide,0.2)
wide.3<-getdiffs(allwide,0.3)
wide.4<-getdiffs(allwide,0.4)
wide.3$diff.psm.4<-abs((wide.3$`0.95`-wide.4$`0.95`)/wide.3$`0.95`)
wide.3$diff.psm.2<-abs((wide.3$`0.95`-wide.2$`0.95`)/wide.3$`0.95`)

#read in the score using chin as metric
chin<-read.csv(here("analysis","results","scores","0.95_0.3_scorechin.csv"))
chinsub<-subset(chin,select=c(ID,score))
chinsub$ID<-as.integer(chinsub$ID)
nsp<-read.csv(here("analysis","results","scores","0.95_0.3_scorensp.csv"))
nspsub<-subset(nsp,select=c(ID,score))
wide.3$ID<-as.integer(wide.3$ID)

chincoho<-left_join(wide.3,chinsub)
chincoho<-chincoho[-which(is.na(chincoho$score)),]
colnames(chincoho)[10]<-"chinscore"
nspcoho<-left_join(wide.3,nspsub)
nspcoho<-nspcoho[-which(is.na(nspcoho$score)),]
colnames(nspcoho)[10]<-"nspscore"
chincoho$diff.chin<-abs((chincoho$`0.95`-chincoho$chinscore)/chincoho$`0.95`)
chincoho$rank<- order(chincoho$`0.95`)
chincoho$rank.chin<- order(chincoho$chinscore)

nspcoho$diff.nsp<-abs((nspcoho$`0.95`-nspcoho$nspscore)/nspcoho$`0.95`)
nspcoho$rank<- order(nspcoho$`0.95`)
nspcoho$rank.nsp<- order(nspcoho$nspscore)

wide.3$`0.95`<-as.numeric(wide.3$`0.95`)
wide.3$rank<- order(wide.3$`0.95`)
wide.3$rank.8<- order(wide.3$`0.8_`)
wide.3$rank.9<- order(wide.3$`0.9_`)

wide.4$rank<- order(wide.4$`0.95`)
wide.2$rank<- order(wide.2$`0.95`)


windows(width=10,height=15)
par(mfrow=c(3,2))
plot(wide.3$rank,wide.3$rank.9)
plot(wide.3$rank,wide.3$rank.8)
plot(wide.3$rank,wide.2$rank)
plot(wide.3$rank,wide.4$rank)
plot(chincoho$rank,chincoho$rank)
plot(nspcoho$rank,nspcoho$rank)
#plot diffs on a single plot
x<-c(1,2,3,4,5,6)
y<-c(mean(wide.3$diff.psm.2),mean(wide.3$diff.psm.4),mean(wide.3$diff.8),mean(wide.4$diff.9),mean(chincoho$diff.chin),mean(nspcoho$diff.nsp))
y.sd<-c(sd(wide.3$diff.psm.2),sd(wide.3$diff.psm.4),sd(wide.3$diff.8),sd(wide.4$diff.9),sd(chincoho$diff.chin),sd(nspcoho$diff.nsp))
windows()
plot(x,y,type="p",pch=16,col="lightblue", ylab="Difference in Prioritization", xlab="",xaxt="n", ylim=c(0,1))
arrows(x,y+y.sd,x,y-y.sd,code=3,length=0,col="gray",lwd=2 )
points(x,y,pch=16,col="lightblue",cex=1.5)
axis(side=1,at=c(1,2,3,4,5,6),labels=c("0.2","0.4","0.8","0.9","Chin","Nsp"))
axis(side=1,at=c(1.5,3.5,5.5),labels=c("psm","alpha","Metric"),line=2, tick=FALSE, lwd=0)

r<-c(cor(wide.3$`0.95`,wide.3$`0.9_`),
     cor(wide.3$`0.95`,wide.3$`0.8_`),
     cor(wide.3$`0.95`,wide.2$`0.95`),
     cor(wide.3$`0.95`,wide.4$`0.95`),
     cor(chincoho$`0.95`,chincoho$chinscore),
     cor(nspcoho$`0.95`,nspcoho$nspscore))
windows()
plot(x,r,type="p",pch=16,col="lightblue", ylab="Correlation of Prioritization Scores", xlab="",xaxt="n", ylim=c(-0.5,1),cex=1.5)
axis(side=1,at=c(1,2,3,4,5,6),labels=c("0.2","0.4","0.8","0.9","Chin","Nsp"))
axis(side=1,at=c(1.5,3.5,5.5),labels=c("psm","alpha","Metric"),line=2, tick=FALSE, lwd=0)

windows(width=10,height=15)
par(mfrow=c(3,2))
plot(wide.3$`0.95`,wide.3$`0.9_`)
plot(wide.3$`0.95`,wide.3$`0.8_`)
plot(wide.3$`0.95`,wide.2$`0.95`)
plot(wide.3$`0.95`,wide.4$`0.95`)
plot(chincoho$`0.95`,chincoho$chinscore)
plot(nspcoho$`0.95`,nspcoho$nspscore)

