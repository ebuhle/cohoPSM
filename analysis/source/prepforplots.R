#preps data files for plotting (in deltaZ_vs_benefits.R or for shiny app in app.R)
## as a mock up of what we may want to do, use mean abundance of spawning salmon as the measure of abundance.
spawnmn<-aggregate(spawn$n,list(spawn$ID,spawn$site),mean)
colnames(spawnmn)<-c("ID","site", "spawn.n")
spawnmn$ID<-as.character(spawnmn$ID)
#get site names associate with ID numbers
siteid<-subset(spatial,select=c(ID,site))
siteid$ID<-as.character(siteid$ID)
spatial_pred$ID<-as.character(spatial_pred$ID)


## merge in spawner data
if(allsites==FALSE){
  psm_pre<-psm_pre[1:51,]
  psm_pre$site<-factor(psm_pre$site)
  d1<-full_join(psm_pre,spawnmn)}
if(allsites==TRUE){
  d1<-full_join(psm_pre,spawnmn)
  spatial<-left_join(spatial_pred,siteid)
  d1<-left_join(d1,siteid)
  #get ID from site column into ID column for sites that don't have it currently.
  d1$ID[52:length(d1$ID)]<-d1$site[52:length(d1$ID)]#this is a problem- there are duplicate ID numbers..
  #for now just remove the rows without a site name 
  d1<-d1[-which(as.numeric(d1$ID)<56 & as.numeric(d1$site)<56),]
  }

#merge the salmonscape data with the spatial data
colnames(salmon)[1]<-"ID"
salmon$ID<-as.character(salmon$ID)
spatial2<-left_join(spatial,salmon)

## add in spatial data
d2<-full_join(d1,spatial2, by="ID")
d<-d2

#plot(psm_pre$Z.mean,psm_pre$p.psm.mean,)
#A cheap way to get Zcrit associated with threshold PSM. 
#Need to make this more robust, use model/include full posterior distribution to get error, etc

Zcrit<-min(d$Z_mean[d$p_psm_mean>input$psm_thresh], na.rm=TRUE)

## Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
## for all sites and select out just the bad sites
d$Zcrit<-Zcrit
d$deltaZ<-d$Zcrit-d$Z_mean


#If desired, use pared down psm_pre file, with just named sites (i.e. sites for which PSM was measured rather than predicted from mdoel)
#if(allsites==FALSE){d<-d[1:51,]}

#add a column for the colors to plot for whether or not site has below threshold psm
d$psmcol<-"darkred" 
d$psmcol[d$p_psm_mean<input$psm_thresh]<-"lightgreen"

