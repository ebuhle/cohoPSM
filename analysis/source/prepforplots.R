#preps data files for plotting (in deltaZ_vs_benefits.R or for shiny app in app.R)
## as a mock up of what we may want to do, use mean abundance of spawning salmon as the measure of abundance.
#spawnmn<-aggregate(spawn$n,list(spawn$ID,spawn$site),mean)
#colnames(spawnmn)<-c("ID","site", "spawn.n")
#spawnmn$ID<-as.character(spawnmn$ID)

#get site names associate with ID numbers
#siteid<-subset(spatial,select=c(ID,site))
#siteid$ID<-as.character(siteid$ID)
#spatial_pred$ID<-as.character(spatial_pred$ID)
#Use only predicted data 

#Sum up the three categories to get total number of meters in which each species is present
salmon$ChinFa_Presence_m<-salmon$Total.FallChin
salmon$ChinSp_Presence_m<-salmon$Total.SprChin
salmon$ChinSu_Presence_m<-salmon$Total.SumChin
salmon$ChumFa_Presence_m<-salmon$Total.FallChum
salmon$ChumWi_Presence_m<-salmon$Total.WinChum
salmon$ChumSu_Presence_m<-salmon$Total.SumChum
salmon$Coho_Presence_m<-salmon$Total.Coho

#Now 
salmon$ChinFa_Presence<-salmon$Total.FallChin
salmon$ChinSp_Presence<-salmon$Total.SprChin
salmon$ChinSu_Presence<-salmon$Total.SumChin
salmon$ChumWi_Presence<-salmon$Total.WinChum
salmon$ChumFa_Presence<-salmon$Total.FallChum
salmon$ChumSu_Presence<-salmon$Total.SumChum
salmon$Coho_Presence<-salmon$Total.Coho

salmon$ChinSp_Presence[!is.na(salmon$ChinSp_Presence) & salmon$ChinSp_Presence>0]<-1
salmon$ChumWi_Presence[!is.na(salmon$ChumWi_Presence) & salmon$ChumWi_Presence>0]<-1
salmon$ChinFa_Presence[!is.na(salmon$ChinFa_Presence) & salmon$ChinFa_Presence>0]<-1
salmon$ChumFa_Presence[!is.na(salmon$ChumFa_Presence) & salmon$ChumFa_Presence>0]<-1
salmon$ChinSu_Presence[!is.na(salmon$ChinSu_Presence) & salmon$ChinSu_Presence>0]<-1
salmon$ChumSu_Presence[!is.na(salmon$ChumSu_Presence) & salmon$ChumSu_Presence>0]<-1
salmon$Coho_Presence[!is.na(salmon$Coho_Presence) & salmon$Coho_Presence>0]<-1
#Add a column for presence of any chum (spr, winter, sum) or chinook (spr, fall, sum)
salmon$Chin_Presence<-NA
salmon$Chin_Presence[salmon$ChinSp_Presence==1]<-1
salmon$Chin_Presence[salmon$ChinSu_Presence==1]<-1
salmon$Chin_Presence[salmon$ChinFa_Presence==1]<-1

#Add a column for presence of any chum (spr, winter, sum) or chinook (spr, fall, sum)
salmon$Chum_Presence<-NA
salmon$Chum_Presence[salmon$ChumWi_Presence==1]<-1
salmon$Chum_Presence[salmon$ChumSu_Presence==1]<-1
salmon$Chum_Presence[salmon$ChumFa_Presence==1]<-1

salmon$nsp_pres<-rowSums(cbind(salmon$Chin_Presence,salmon$Coho_Presence,salmon$Chum_Presence),na.rm=TRUE)

## merge spatial, psm, and salmon information
##do not merge in spawner data to predicted stuff
if(predsites==FALSE){
  psm_pre<-psm_pre[1:51,]
  psm_pre$site<-factor(psm_pre$site)
  d3<-full_join(psm_pre,spawnmn)}
if(predsites==TRUE){
  goo<-psm_pre[52:dim(psm_pre)[1],]
  colnames(goo)[1]<-"ID"
  goo$ID<-as.integer(goo$ID)
  d1<-left_join(goo,spatial_pred)
    }

#merge the salmonscape data with the spatial data
colnames(salmon)[1]<-"ID"
salmon$ID<-as.integer(salmon$ID)

## add in spatial data
d<-full_join(d1,salmon, by="ID")

#plot(psm_pre$Z.mean,psm_pre$p.psm.mean,)
#A cheap way to get Zcrit associated with threshold PSM. 
#Need to make this more robust, use model/include full posterior distribution to get error, etc

