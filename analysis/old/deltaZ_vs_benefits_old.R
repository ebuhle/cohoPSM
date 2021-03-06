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
#######################################
###Below is old code that:
###1) uses delta z instead of Z to plot the points
###2) divides into restoration and conservation- we decided not to do this
#######################################

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
plot(d$Z_mean, d$spawn.n, cex=2,cex.lab=1.2,cex.axis=1.2,xlab=expression(paste("Change in urbanization required (",Delta,"Z)", sep="")), ylab= "Benefit = spawner abundance", type="p", pch=19, col=d$psmcol)
## Step 4. Compare to other benefits of interest- abundance or presence of other salmon species for all sites, human pops, stream biodiversity. Can blake get these?
text(d$Z_mean, d$spawn.n, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(0,0,min(d$deltaZ, na.rm=TRUE),min(d$deltaZ, na.rm=TRUE)),
        c(0,max(d$spawn.n, na.rm=TRUE),max(d$spawn.n, na.rm=TRUE),0),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
## Step 4. Compare to other benefits of interest- abundance or presence of other salmon species for all sites, human pops, stream biodiversity. Can blake get these?
text(d$Z_mean, d$spawn.n, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(0,0,Zcrit,Zcrit),
        c(0,max(d$spawn.n, na.rm=TRUE),max(d$spawn.n, na.rm=TRUE),0),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)
abline(v=Zcrit, lty=2, lwd=2, col="blue")

mtext(side=1,"high",line=4,adj=1)
mtext(side=1,"low",line=4,adj=0)

head(d)

#plot change in Z vs benefit of number of people living nearby (don't have this data right now, so making it up as a function of Z with some noist)
#plot(psm_pre3$Z.mean,pop.pred)#check
plot(d$deltaZ,d$pop_census, cex=2,cex.lab=1.2,cex.axis=1.2,xlab=expression(paste("Change in urbanization required (",Delta,"Z)", sep="")), ylab= "Benefit = number of people nearby", type="p", pch=19, col=d$psmcol)
text(d$deltaZ, d$pop_census, labels=as.numeric(as.factor(d$site)),cex=0.8, font=2)
polygon(c(0,0,min(d$deltaZ, na.rm=TRUE),min(d$deltaZ, na.rm=TRUE)),
        c(0,max(d$pop_census, na.rm=TRUE),max(d$pop_census),0),col=adjustcolor("salmon",alpha.f=0.5),
        border=NA)

dev.off()

#plot change in Z vs benefit (first  benefit=spawner abundance)
dev.new(height=5,width=10)
par(mfrow=c(1,2))
#restoration sites
#sites with greater than threshold psm- in need of restoration. 
r<-d[d$p_psm_mean>=input$psm_thresh,]
#standardizgin the effort and the benefit, so that they are equally weighted...we can decide if we want to weight things differently.

r$deltaZ.stan<-(r$deltaZ-mean(r$deltaZ))/sd(r$deltaZ)
r$spawn.n.stan<-(r$spawn.n-mean(r$spawn.n))/sd(r$spawn.n)
rxy<-subset(r,select=c(deltaZ.stan,spawn.n.stan))
rest.score<-as.matrix(dist(rbind(c(0,max(rxy$deltaZ.stan)),rxy), method="euclidean"))[1,-1]
rxy<-cbind(r$site,r$deltaZ,r$spawn.n,rxy,rest.score)
rxy<-rxy[order(rxy$rest.score),]
myPalette <- colorRampPalette(brewer.pal(9, "RdYlGn")) #### Gives us a heat map look
cols = myPalette(length(rest.score))
rxy<- data.frame(cbind(rxy,cols))
colnames(rxy)[1:3]<-c("site","deltaZ","benefit")
plot(rxy$deltaZ, rxy$benefit, cex=2,cex.lab=1.2,cex.axis=1.2,xlim=range(rxy$deltaZ),xlab="Effort", ylab= "Benefit = spawner abundance", type="p", pch=19, col=as.character(rxy$cols))
mtext(side=1,expression(paste("Change in urbanization required (",Delta,"Z)", sep="")),line=4)
mtext(side=3,"Goal=Restoration",line=0, cex=1.2)
mtext(side=1,"high",line=4,adj=0)
mtext(side=1,"low",line=4,adj=1)
#mtext(side=1,"Highest priority", line=-20, adj=.9, cex=1.2, col="darkgreen")
legend("topright", legend=c("Highest priority","Lowest priority"), pch=19,col=c(cols[length(cols)],cols[1]))
text(rxy$deltaZ, rxy$benefit, labels=as.numeric(as.factor(rxy$site)),cex=0.8, font=2)

spawnscore<-as.data.frame(cbind(rxy$site,rxy$rest.score))
colnames(spawnscore)<-c("site","restscore.spawners")

#now try for a different benefit: # of people nearby
r$pop_census.stan<-(r$pop_census-mean(r$pop_lscan))/sd(r$pop_census)
rxy<-subset(r,select=c(deltaZ.stan,pop_census.stan))
rest.score<-as.matrix(dist(rbind(c(0,max(rxy$deltaZ.stan)),rxy), method="euclidean"))[1,-1]
rxy<-cbind(r$site,r$deltaZ,r$pop_census,rxy,rest.score)
rxy<-rxy[order(rxy$rest.score),]
myPalette <- colorRampPalette(brewer.pal(9, "RdYlGn")) #### Gives us a heat map look
cols = myPalette(length(rest.score))
rxy<- data.frame(cbind(rxy,cols))
colnames(rxy)[1:3]<-c("site","deltaZ","benefit")
plot(rxy$deltaZ, rxy$benefit, cex=2,cex.lab=1.2,cex.axis=1.2,xlim=range(rxy$deltaZ),xlab="Effort", ylab= "Benefit = # of People Living Nearby", type="p", pch=19, col=as.character(rxy$cols))
mtext(side=1,expression(paste("Change in urbanization required (",Delta,"Z)", sep="")),line=4)
mtext(side=3,"Goal=Restoration",line=0, cex=1.2)
mtext(side=1,"high",line=4,adj=0)
mtext(side=1,"low",line=4,adj=1)
text(rxy$deltaZ, rxy$benefit, labels=as.numeric(as.factor(rxy$site)),cex=0.8, font=2)

peoplescore<-as.data.frame(cbind(rxy$site,rxy$rest.score))
colnames(peoplescore)<-c("site","restscore.people")
rest.scores<-full_join(spawnscore,peoplescore)
write.csv(rest.scores,here("analysis","output","restscores.csv"), row.names=FALSE)
#conservation sites= those with below threshold psm
dev.new(height=5,width=10)
par(mfrow=c(1,2))
c<-d[d$p_psm_mean<input$psm_thresh,]
c$deltaZ.stan<-(c$deltaZ-mean(c$deltaZ))/sd(c$deltaZ)
c$spawn.n.stan<-(c$spawn.n-mean(c$spawn.n))/sd(c$spawn.n)
cxy<-subset(c,select=c(deltaZ.stan,spawn.n.stan))
con.score<-as.matrix(dist(rbind(c(0,min(cxy$deltaZ.stan)),cxy), method="euclidean"))[1,-1]
cxy<-cbind(c$site,c$deltaZ,c$spawn.n,cxy,con.score)
cxy<-cxy[order(cxy$con.score),]
myPalette <- colorRampPalette(brewer.pal(9, "RdYlGn")) #### Gives us a heat map look
cols = myPalette(length(con.score))
cxy<- data.frame(cbind(cxy,cols))
colnames(cxy)[1:3]<-c("site","deltaZ","spawn.n")
plot(cxy$deltaZ, cxy$spawn.n, cex=2,cex.lab=1.2,cex.axis=1.2,xlim=range(cxy$deltaZ),xlab="?", ylab= "Benefit = spawner abundance", type="p", pch=19, col=as.character(cxy$cols))
mtext(side=1,expression(paste("Change in urbanization required (",Delta,"Z)", sep="")),line=4)
mtext(side=3,"Goal=Conservation",line=0, cex=1.2)
#mtext(side=1,"high",line=4,adj=0)
#mtext(side=1,"low",line=4,adj=1)
text(cxy$deltaZ, cxy$spawn.n, labels=as.numeric(as.factor(cxy$site)),cex=0.8, font=2)
spawnscore.cons<-as.data.frame(cbind(cxy$site,cxy$con.score))
colnames(spawnscore.cons)<-c("site","conscore.spawners")
#now add human pop as a benefit
c$pop_census.stan<-(c$pop_census-mean(c$pop_census))/sd(c$pop_census)
cxy<-subset(c,select=c(deltaZ.stan,pop_census.stan))
con.score<-as.matrix(dist(rbind(c(0,min(cxy$deltaZ.stan)),cxy), method="euclidean"))[1,-1]
cxy<-cbind(c$site,c$deltaZ,c$pop_census,cxy,con.score)
cxy<-cxy[order(cxy$con.score),]
myPalette <- colorRampPalette(brewer.pal(9, "RdYlGn")) #### Gives us a heat map look
cols = myPalette(length(con.score))
cxy<- data.frame(cbind(cxy,cols))
colnames(cxy)[1:3]<-c("site","deltaZ","pop_census")
plot(cxy$deltaZ, cxy$pop_census, cex=2,cex.lab=1.2,cex.axis=1.2,xlim=range(cxy$deltaZ),xlab="?", ylab= "Benefit = # humans living nearby", type="p", pch=19, col=as.character(cxy$cols))
mtext(side=1,expression(paste("Change in urbanization required (",Delta,"Z)", sep="")),line=4)
mtext(side=3,"Goal=Conservation",line=0, cex=1.2)
#mtext(side=1,"high",line=4,adj=0)
#mtext(side=1,"low",line=4,adj=1)
text(cxy$deltaZ, cxy$pop_census, labels=as.numeric(as.factor(cxy$site)),cex=0.8, font=2)
legend("topright", legend=c("Highest priority","Lowest priority"), pch=19,col=c(cols[length(cols)],cols[1]))




## I looked into salmon abundance data from DWF stock inventory pops- no "bad sites" included. some good sites included...look into this more
badsites<-psm_pre[psm_pre$Z_mean>psm_thresh,]
sitenames<-unique(psm_pre$site)[1:51]
badsites<-badsites[1:24,]#just the sites with names

abund<-read.csv(here("data","WDFW-Salmonid_Stock_Inventory_Populations.csv"), header=TRUE)
sort(unique(abund$Population.Name))
dim(badsites)
tail(badsites)
badsites$site
#Questions:
#What are the numbered sites?