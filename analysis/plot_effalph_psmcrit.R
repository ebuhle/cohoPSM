  #make a table andplot showing how change sin alpha and psmcrit affect prioritization
  #started by ailene
  #dec 22, 2019
  
  ## housekeeping
  rm(list=ls()) 
options(device = ifelse(.Platform$OS.type == "windows", "windows", "quartz"))
#options(stringsAsFactors = FALSE)
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
files.95.3<-files[grep("0.95_0.3_",files)]
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
head(alphpsmcrit)
#convert allfiles from long format to wide format
allfiles<-allfiles[,-which(colnames(allfiles)=="cols")]
allfiles$ID<-factor(allfiles$ID)
allfiles$alpha<-factor(allfiles$alpha)

allfiles$psmcrit<-factor(allfiles$psmcrit)
allfiles<-left_join(allfiles,allconsrest)

allfiles<-allfiles[order(allfiles$ID,allfiles$alpha,allfiles$psmcrit),]
allfiles$action<-"cons"
allfiles$action[allfiles$delta_Z<0]<-"rest"
allsub<-subset(allfiles,select=c(ID,priority,alpha,psmcrit))
allwide<-spread(allsub,alpha,priority)
allwide$psmcrit<-as.numeric(as.character(allwide$psmcrit))

wide.1<-allwide[allwide$psmcrit==0.1,]
wide.2<-allwide[allwide$psmcrit==0.2,]
wide.4<-allwide[allwide$psmcrit==0.4,]
wide.5<-allwide[allwide$psmcrit==0.5,]
wide.15<-allwide[allwide$psmcrit==0.15,]
wide.25<-allwide[allwide$psmcrit==0.25,]
wide.35<-allwide[allwide$psmcrit==0.35,]
wide.45<-allwide[allwide$psmcrit==0.45,]
wide.55<-allwide[allwide$psmcrit==0.55,]

getdiffs<-function(allwide,psmcrit){
  w<-allwide[allwide$psmcrit==psmcrit,]
  w<-w[order(w$ID),]
  w$diff.99<-abs(w$`0.95`-w$`0.99`)
  w$diff.95<-abs(w$`0.95`-w$`0.95`)
  w$diff.9<-abs(w$`0.95`-w$`0.9_`)
  w$diff.85<-abs(w$`0.95`-w$`0.85`)
  w$diff.8<-abs(w$`0.95`-w$`0.8_`)
  w$diff.75<-abs(w$`0.95`-w$`0.75`)
  w$diff.7<-abs(w$`0.95`-w$`0.7_`)
  
  return(w) 
}

wide.3<-getdiffs(allwide,.3)
wide.3$diff.psm.1<-abs(wide.3$`0.95`-wide.1$`0.95`)
wide.3$diff.psm.15<-abs(wide.3$`0.95`-wide.15$`0.95`)
wide.3$diff.psm.2<-abs(wide.3$`0.95`-wide.2$`0.95`)
wide.3$diff.psm.25<-abs(wide.3$`0.95`-wide.25$`0.95`)
wide.3$diff.psm.35<-abs(wide.3$`0.95`-wide.35$`0.95`)
wide.3$diff.psm.4<-abs(wide.3$`0.95`-wide.4$`0.95`)
wide.3$diff.psm.45<-abs(wide.3$`0.95`-wide.45$`0.95`)
wide.3$diff.psm.5<-abs(wide.3$`0.95`-wide.5$`0.95`)

#read in the score using chin as metric
chin<-read.csv(here("analysis","results","scores","0.95_0.3_scorechin.csv"))
chinsub<-subset(chin,select=c(ID,priority))
chinsub$ID<-as.integer(chinsub$ID)
nsp<-read.csv(here("analysis","results","scores","0.95_0.3_scorensp.csv"))
nspsub<-subset(nsp,select=c(ID,priority))
wide.3$ID<-as.integer(wide.3$ID)

chincoho<-left_join(wide.3,chinsub)
chincoho<-chincoho[-which(is.na(chincoho$priority)),]
colnames(chincoho)[which(colnames(chincoho)=="priority")]<-"chinscore"
nspcoho<-left_join(wide.3,nspsub)
nspcoho<-nspcoho[-which(is.na(nspcoho$priority)),]
colnames(nspcoho)[which(colnames(nspcoho)=="priority")]<-"nspscore"
chincoho$diff.chin<-abs(chincoho$`0.95`-chincoho$chinscore)

chincoho$rank<- order(chincoho$`0.95`)
chincoho$rank.chin<- order(chincoho$chinscore)

nspcoho$diff.nsp<-abs(nspcoho$`0.95`-nspcoho$nspscore)
nspcoho$rank<- order(nspcoho$`0.95`)
nspcoho$rank.nsp<- order(nspcoho$nspscore)

wide.3$`0.95`<-as.numeric(wide.3$`0.95`)
wide.3$rank<- order(wide.3$`0.95`)
wide.3$rank.8<- order(wide.3$`0.8_`)
wide.3$rank.9<- order(wide.3$`0.9_`)

wide.4$rank<- order(wide.4$`0.95`)
wide.2$rank<- order(wide.2$`0.95`)
#plot diffs on 3 panels
psm_crit_vals <- seq(0.1, 0.5, 0.05)
alpha_vals <- seq(0.7, 0.99, 0.05)

png(here("analysis","results","figures","metricuncertaintycomp.png"),height = 500,width = 1000)

par(mfrow=c(1,2))
x<-psm_crit_vals
y<-c(mean(wide.3$diff.psm.1),
     mean(wide.3$diff.psm.15),
     mean(wide.3$diff.psm.2),
     mean(wide.3$diff.psm.25),
     mean(wide.3$diff.psm.3),
     mean(wide.3$diff.psm.35),
     mean(wide.3$diff.psm.4),
     mean(wide.3$diff.psm.45),
     mean(wide.3$diff.psm.5))
y.sd<-c(sd(wide.3$diff.psm.1),
        sd(wide.3$diff.psm.15),
        sd(wide.3$diff.psm.2),
        sd(wide.3$diff.psm.25),
        sd(wide.3$diff.psm.3),
        sd(wide.3$diff.psm.35),
        sd(wide.3$diff.psm.4),
        sd(wide.3$diff.psm.45),
        sd(wide.3$diff.psm.5))

plot(x,y,type="p",pch=16,col="darkgray", ylab="Change in prioritization rank", xlab=expression(paste("PSM"[crit])), ylim=c(0,450), bty="l",cex.lab= 1.5,cex.axis = 1.5 )
arrows(x,y+y.sd,x,y-y.sd,code=3,length=0,col="lightgray",lwd=2 )
points(x,y,pch=16,col="darkgray",cex=2)
points(0.3, mean(chincoho$diff.chin), pch = 16,cex = 2,col = "darkblue" )
arrows(0.3,mean(chincoho$diff.chin)+sd(chincoho$diff.chin),0.3,mean(chincoho$diff.chin)-sd(chincoho$diff.chin),code=3,length=0,col="darkblue",lwd=2 )
mtext(side= 3,"A)", line = 1, adj=0,cex=2)
#alpha/uncertainty
x=c(0.70,0.75,0.80,0.85,0.90,0.95,0.99)
y<-c(mean(wide.3$diff.7),
       mean(wide.3$diff.75),
       mean(wide.3$diff.8),
       mean(wide.3$diff.85),
       mean(wide.3$diff.9),
       mean(wide.3$diff.95),
       mean(wide.3$diff.99))
y.sd<-c(sd(wide.3$diff.7),
     sd(wide.3$diff.75),
     sd(wide.3$diff.8),
     sd(wide.3$diff.85),
     sd(wide.3$diff.9),
     sd(wide.3$diff.95),
     sd(wide.3$diff.99))
plot(x,y,type="p",pch=16,col="darkgray", ylab="Change in Prioritization Rank", xlab="Uncertainty",ylim=c(0,450), bty="l",cex.lab= 1.5,cex.axis = 1.5)
arrows(x,y+y.sd,x,y-y.sd,code=3,length=0,col="gray",lwd=2 )
points(x,y,pch=16,col="darkgray",cex=2)
points(0.95, mean(chincoho$diff.chin), pch = 16,cex = 2,col = "darkblue" )
arrows(0.95,mean(chincoho$diff.chin)+sd(chincoho$diff.chin),0.95,mean(chincoho$diff.chin)-sd(chincoho$diff.chin),code=3,length=0,col="darkblue",lwd=2 )
mtext(side= 3,"B)", line = 1, adj=0,cex=2)

#metric
# x=c(.5,1.5)
# y<-c(mean(chincoho$diff.chin),
#      mean(nspcoho$diff.nsp))
# y.sd<-c(sd(chincoho$diff.chin),
#         sd(nspcoho$diff.nsp))
# plot(x,y,type="p",pch=16,col="darkgray", ylab="Difference in Prioritization (Proportion)", xlab="Metric", xaxt="n",xlim=c(-1,3),ylim=c(0,400))
# arrows(x,y+y.sd,x,y-y.sd,code=3,length=0,col="gray",lwd=2)
# points(x,y,pch=16,col="darkgray",cex=2)
# axis(side=1,at=c(0.5, 1.5),labels=c("Chinook","Nsp"),line=0,tick=FALSE, lwd=0)
# 
# 
dev.off()
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

