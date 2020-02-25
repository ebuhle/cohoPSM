#Make Table S1 For coho supplement 
#Table S1 includes the scores and delta zs and rankings of each subbasin in our study
library(here)
library(dplyr)
#read in score files
coho<-read.csv(here("analysis","results","scores","0.95_0.3_score.csv"))
chin<-read.csv(here("analysis","results","scores","0.95_0.3_scorechin.csv"))
nsp<-read.csv(here("analysis","results","scores","0.95_0.3_scorensp.csv"))
coho$conservation.action<-"preservation"
coho$conservation.action[coho$delta_Z<0]<-"restoration"
#get number of sites in each
table(coho$conservation.action)

coho<-coho[order(coho$conservation.action,coho$score),]
coho$coho.rank[coho$conservation.action=="preservation"]<-seq(1:length(coho$conservation.action[coho$conservation.action=="preservation"]))
coho$coho.rank[coho$conservation.action=="restoration"]<-seq(1:length(coho$conservation.action[coho$conservation.action=="restoration"]))
coho<-subset(coho,select=c(ID,conservation.action,benefit,coho.rank))
colnames(coho)[3]<-"Coho_Habitat_m"

chin$conservation.action<-"preservation"
chin$conservation.action[chin$delta_Z<0]<-"restoration"
chin<-chin[order(chin$conservation.action,chin$score),]
chin$chin.rank[chin$conservation.action=="preservation"]<-seq(1:length(chin$conservation.action[chin$conservation.action=="preservation"]))
chin$chin.rank[chin$conservation.action=="restoration"]<-seq(1:length(chin$conservation.action[chin$conservation.action=="restoration"]))
table(chin$conservation.action)
chin<-subset(chin,select=c(ID,conservation.action,benefit,chin.rank))
colnames(chin)[3]<-"Chin_Habitat_m"

nsp$conservation.action<-"preservation"
nsp$conservation.action[nsp$delta_Z<0]<-"restoration"
nsp<-nsp[order(nsp$conservation.action,nsp$score),]
nsp$nsp.rank[nsp$conservation.action=="preservation"]<-seq(1:length(nsp$conservation.action[nsp$conservation.action=="preservation"]))
nsp$nsp.rank[nsp$conservation.action=="restoration"]<-seq(1:length(nsp$conservation.action[nsp$conservation.action=="restoration"]))
table(nsp$conservation.action)
nsp<-subset(nsp,select=c(ID,conservation.action,benefit,nsp.rank))
colnames(nsp)[3]<-"nsp"
head(nsp)
head(chin)
head(coho)
tab<-left_join(coho,chin)
tab2<-left_join(tab,nsp)
tab2<-tab2[order(tab2$ID),]
write.csv(tab2,"results/prioritization_TableS1.csv", row.names=FALSE)

