## libraries
library(dplyr)
library(shiny)
library(RColorBrewer)
## Read in the data/model estimates and join
psm_pre<-read.table("../analysis/results/PSM_predictions.txt",header=TRUE)
spawn<-read.csv("../data/spawner_data.csv", header=TRUE)
spatial<-read.csv("../data/spatial_data.csv", header=TRUE)

spawnmn<-aggregate(spawn$n,list(spawn$site),mean)
colnames(spawnmn)<-c("site", "spawn.n")
spatial$site<-as.factor(spatial$site)
allsites=FALSE

allsites==FALSE
#eventually, use spatial_data_predict.csv instead of spatial_data.csv

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

#eventually use salmon scape data to add salmon species:
d$salmonsp.n<-1

#add other benefits?

dsubs<-subset(d,select=c(site,Z.mean,p.psm.mean,ID,watershed,spawn.n,salmonsp.n,pop_census)) 
colnames(dsubs)[6:8]<-c("Num.coho.spawners","Num.salmon.spp","Num.people.nearby")

ui <- pageWithSidebar(
  
  # App title ----
 headerPanel("Prioritizing Coho Sites"),
  
  # Sidebar layout with input and output definitions ----
  #sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
   
      
      # Input: Selector for variable to plot against mpg ----
      selectInput("goal", "Goal:", 
                  c("Conservation","Restoration"),
                  selected="Restoration"),
      
      # Input: Selector for variable to plot against mpg ----
      selectInput("benefit", "Benefit:", 
                  names(dsubs[,6:8]),
                  selected=names(dsubs[,6:8])[1]),
      
      # Input: Slider for selecting critical PSM threshold ----
      sliderInput(inputId = "psm_thresh",
                  label = "Pre-Spawn Mortality (PSM) Threshold",
                  min = 0,
                  max = 1.0,
                  value = 0.25)
      ),
      
   
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plot PSM vs. Z ----
      plotOutput(outputId = "plot_psm_z"),
      
      # Output: Plot selected benefit vs deltaZ
      plotOutput(outputId = "plot_benefit")
    )
  )


# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  # Plot PSM vs. Urbanization (z) ----
  # with requested Threshold PSM highlighted
  # This generates plot that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$plot_psm_z <- renderPlot({
    
    #A cheap way to get Zcrit associated with threshold PSM. 
    #Need to make this more robust, use model/include full posterior distribution to get error, etc
    
    Zcrit<-min(dsubs$Z.mean[dsubs$p.psm.mean>input$psm_thresh])
    
    ## Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    ## for all sites and select out just the bad sites
    dsubs$Zcrit<-Zcrit
    dsubs$deltaZ<-dsubs$Zcrit-d$Z.mean
    
    #add a column for the colors to plot for whether or not site has below threshold psm
    dsubs$psmcol<-"darkred" 
    dsubs$psmcol[dsubs$p.psm.mean<input$psm_thresh]<-"lightgreen"
   
    plot(dsubs$p.psm.mean, dsubs$Z.mean, pch=19,col=dsubs$psmcol, cex.lab=1.2,cex.axis=1.2,cex=2, ylab="Urbanization effects (Z)", xlab= "Mean Pre-Spawn Mortality")
    
    abline(v=input$psm_thresh, lty=2, lwd=2)
    
    text(input$psm_thresh+.02,min(d$Z.mean),label="PSM threshold", cex=1.2)
    
    abline(h=Zcrit, lty=2, lwd=2, col="blue")
    
    polygon(c(input$psm_thresh,1,1,input$psm_thresh),c(Zcrit,Zcrit,max(d$Z.mean)+.5,max(d$Z.mean)+.5),
            col=adjustcolor("salmon",alpha.f=0.5),
            border=NA)
    
    text(.02,Zcrit+.04,label="Zcrit", col="blue",cex=1.2)
    
    
  })
  
  # selectedData <- reactive({
  #   dsubs[, c("site",input$benefit,"deltaZ","p.psm.mean")]
  # })
  output$plot_benefit<- renderPlot({
    dsubs2<-dsubs[, c("site",input$benefit,"deltaZ","p.psm.mean")]
    r<-dsubs[dsubs2$p.psm.mean>=input$psm_thresh,]
    #standardizing the effort and the benefit, so that they are equally weighted...we can decide if we want to weight things differently.
    r$deltaZ.stan<-(r$deltaZ-mean(r$deltaZ))/sd(r$deltaZ)
   
    r$benefit.stan<-(r[,2]-mean(r[,2]))/sd(r[,2])
    rxy<-r[,5:6]
    rest.score<-as.matrix(dist(rbind(c(0,max(rxy$deltaZ.stan)),rxy), method="euclidean"))[1,-1]
    rxy<-cbind(r[, c("site","deltaZ",input$benefit)],rxy,rest.score)
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
  })
  
}

shinyApp(ui = ui, server = server)
