#preps data files for plotting (in deltaZ_vs_benefits.R or for shiny app in app.R)
## as a mock up of what we may want to do, use mean abundance of spawning salmon as the measure of abundance.
spawnmn<-aggregate(spawn$n,list(spawn$site),mean)
colnames(spawnmn)<-c("site", "spawn.n")
spatial$site<-as.factor(spatial$site)
## merge in spawner data
if(allsites==FALSE){
  psm_pre<-psm_pre[1:51,]
  psm_pre$site<-factor(psm_pre$site)
  d<-full_join(psm_pre,spawnmn)}
if(allsites==TRUE){
  d<-psm_pre
  d<-full_join(psm_pre,spawnmn)}

## add in spatial data
d<-full_join(d,spatial)


#plot(psm_pre$Z.mean,psm_pre$p.psm.mean,)
#A cheap way to get Zcrit associated with threshold PSM. 
#Need to make this more robust, use model/include full posterior distribution to get error, etc

Zcrit<-min(d$Z_mean[d$p_psm_mean>input$psm_thresh])

## Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
## for all sites and select out just the bad sites
d$Zcrit<-Zcrit
d$deltaZ<-d$Zcrit-d$Z_mean


#If desired, use pared down psm_pre file, with just named sites (i.e. sites for which PSM was measured rather than predicted from mdoel)
#if(allsites==FALSE){d<-d[1:51,]}

#add a column for the colors to plot for whether or not site has below threshold psm
d$psmcol<-"darkred" 
d$psmcol[d$p_psm_mean<input$psm_thresh]<-"lightgreen"

