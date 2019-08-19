#==================================================================
# SETUP
#==================================================================

library(dplyr)
library(yarrr)
library(vioplot)
library(colourvalues)
library(viridis)
library(shape)
library(rstan)
library(Hmisc)
library(matrixStats)
library(here)
library(leaflet)
library(rgdal)

# load functions
source(here("analysis","functions","stan_mean.R"))
source(here("analysis","functions","extract1.R"))
source(here("analysis","functions","sem_psm_predict.R"))
source(here("analysis","functions","sem_lulc_predict.R"))
source(here("analysis","functions","sem_z_crit.R"))
source(here("analysis","functions","vioplot2.R"))

# read and wrangle data
source(here("analysis","cohoPSM1_data.R"))  

# read geospatial data/model estimates (use only predicted attributes for now)
map_pred <- readOGR(here("data","geospatial","predictions_with_all_GIS_data_Albers.shp"))

# load previously saved stanfit objects
if(file.exists(here("analysis","results","stan_psm_all.RData"))) # full SEM incl WADOE basins
  load(here("analysis","results","stan_psm_all.RData"))

# Step 2: Choices: select the threshold psm and you want to use, 
# and select all IDs or only IDs for which we have PSM data (rather than predicted PSM)
input <- as.data.frame(NA)
input$psm_thresh <- 0.25
input$attribute<-"Coho_Presence_m"

predsites <- TRUE #if false, selects out only IDs with PSM calculated from field data, rather than IDs with predicted PSM

## Step 3: combine all the data and prep for plotting calculate mean spawner abundance by ID, across years
## combined data file with things we want to plot is called "d"
attrs <- c("coho_presence_m","N_spp","fallchin_presence_m")

ui <- pageWithSidebar(
  # App title
  headerPanel("Prioritizing Coho Sites"),
  
  # Sidebar panel for inputs
  sidebarPanel(
    # Input: slider for critical PSM threshold
    sliderInput(inputId = "psm_crit",
                label = "PSM threshold",
                min = 0,
                max = 1.0,
                value = 0.3),
    
    # Input: slider for confidence level
    sliderInput(inputId = "alpha",
                label = "Confidence level",
                min = 0,
                max = 1,
                value = 0.9),
    
    # Input: selector for attribute to plot against z
    selectInput(inputId = "attribute", 
                label = "Attribute", 
                attrs,
                selected = attrs[1]),
    
    # Input: selector for focal site
    selectInput(inputId = "show_site",
                label = "Highlight site",
                levels(spatial_data_pre$site),
                selected = levels(spatial_data_pre$site)[1])
  ),
  
  
  
  # Main panel for displaying outputs ----
  mainPanel(
    
    # Output: Plot PSM vs. Z ----
    plotOutput(outputId = "plot_psm_z"),
    
    # Output: Plot selected attribute vs Z
    plotOutput(outputId = "plot_attribute"),
    #output: map of sites
    leafletOutput("mymap",height=1000)
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
  
  
  # Plot PSM vs. Urbanization (z) ----
  # with requested Threshold PSM highlighted
  # This generates plot that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # Return the requested dataset ----
  selectedData <- reactive({
    Zcrit<-min(d$Z_mean[d$p_psm_mean>input$psm_thresh], na.rm=TRUE)
    #Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    # for all sites and select out just the bad sites
    #add a column for the colors to plot for whether or not site has below threshold psm
    attribute<-d[,which(colnames(d)==input$attribute)]
    d$attribut_stan<-(attribute-mean(attribute,na.rm=TRUE))/sd(attribute,na.rm=TRUE)
    dxy<-subset(d,select=c(Z_mean,attribut_stan))
    score<-as.matrix(dist(rbind(c(Zcrit,max(dxy$attribut_stan,na.rm=TRUE)),dxy), method="euclidean"))[1,-1]
    dxy<-cbind(d$ID,dxy,d[,which(colnames(d)==input$attribute)],score)
    dxy<-dxy[-which(is.na(dxy$attribut_stan)),]
    dxy<-dxy[order(dxy$score),]
    myPalette <- colorRampPalette(brewer.pal(9, "RdYlBu")) #### Gives us a heat map look
    cols = rev(myPalette(length(dxy$score)))
    dxy<- data.frame(cbind(dxy,cols))
    colnames(dxy)[1:4]<-c("ID","Z","benefit.stan","benefit")
  })
  # 2. Its output type is a plot
  output$plot_psm_z <- renderPlot({
    
    #A cheap way to get Zcrit associated with threshold PSM. 
    #Need to make this more robust, use model/include full posterior distribution to get error, etc
    
    Zcrit<-min(d$Z_mean[d$p_psm_mean>input$psm_thresh], na.rm=TRUE)
    #Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    # for all sites and select out just the bad sites
    d$Zcrit<-Zcrit
    d$deltaZ<-d$Zcrit-d$Z_mean
    #add a column for the colors to plot for whether or not site has below threshold psm
    d$psmshape<-24 
    d$psmshape[d$p_psm_mean<input$psm_thresh]<-19
    d$psmcol<-"black"
    d$psmcol[d$p_psm_mean<input$psm_thresh]<-"gray"
    d$psmbg<-"gray"
    d$psmbg[d$p_psm_mean<input$psm_thresh]<-NA
    
    ## Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    ## for all sites and select out just the bad sites
    #plot(dsubs$p_psm_mean, dsubs$Z_mean, pch=d$shape,cex.lab=1.2,cex.axis=1.2,cex=2, ylab="Urbanization effects (Z)", xlab= "Mean Pre-Spawn Mortality")
    plot(d$Z_mean,d$p_psm_mean,pch=d$psmshape, col=d$psmcol, bg=d$psmbg,cex.lab=1.2,cex.axis=1.2,cex=1.52, ylim=c(0,1),
         xlab="Urbanization score (Z)", ylab= "Pre-Spawn Mortality")
    abline(h=input$psm_thresh, lty=2, lwd=2)
    text(min(d$Z_mean, na.rm=TRUE)+1,input$psm_thresh+.02,label="PSM threshold", cex=1.2)
    abline(v=Zcrit, lty=2, lwd=2, col="blue")
    
    polygon(c(input$psm_thresh,1,1,input$psm_thresh),c(Zcrit,Zcrit,max(d$Z.mean)+.5,max(d$Z.mean)+.5),
            col=adjustcolor("salmon",alpha.f=0.5),
            border=NA)
    
    text(.02,Zcrit+.04,label="Zcrit", col="blue",cex=1.2)
    legend("topleft", legend=c("Above threshold","Below threshold"), pch=c(24,19),bg=c("gray",NA),col=c("black","gray"), cex=.8, bty="n")
    
  })
  
  output$plot_attribute <- renderPlot({
    
    #A cheap way to get Zcrit associated with threshold PSM. 
    #Need to make this more robust, use model/include full posterior distribution to get error, etc
    
    Zcrit<-min(d$Z_mean[d$p_psm_mean>input$psm_thresh], na.rm=TRUE)
    #Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    # for all sites and select out just the bad sites
    d$Zcrit<-Zcrit
    d$deltaZ<-d$Zcrit-d$Z_mean
    #add a column for the colors to plot for whether or not site has below threshold psm
    d$psmshape<-17 
    d$psmshape[d$p_psm_mean<input$psm_thresh]<-19
    Zcrit<-min(d$Z_mean[d$p_psm_mean>input$psm_thresh], na.rm=TRUE)
    #Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    # for all sites and select out just the bad sites
    d$Zcrit<-Zcrit
    d$deltaZ<-d$Zcrit-d$Z_mean
    #add a column for the colors to plot for whether or not site has below threshold psm
    d$psmshape<-17 
    d$psmshape[d$p_psm_mean<input$psm_thresh]<-19
    attribute<-d[,which(colnames(d)==input$attribute)]
    d$attribut_stan<-(attribute-mean(attribute,na.rm=TRUE))/sd(attribute,na.rm=TRUE)
    dxy<-subset(d,select=c(Z_mean,attribut_stan))
    score<-as.matrix(dist(rbind(c(Zcrit,max(dxy$attribut_stan,na.rm=TRUE)),dxy), method="euclidean"))[1,-1]
    dxy<-cbind(d$ID,dxy,d[,which(colnames(d)==input$attribute)],score)
    dxy<-dxy[-which(is.na(dxy$attribut_stan)),]
    dxy<-dxy[order(dxy$score),]
    myPalette <- colorRampPalette(brewer.pal(9, "RdYlBu")) #### Gives us a heat map look
    cols = rev(myPalette(length(dxy$score)))
    dxy<- data.frame(cbind(dxy,cols))
    colnames(dxy)[1:4]<-c("ID","Z","benefit.stan","benefit")
    
    plot(dxy$Z,dxy$benefit, cex=1.5,cex.lab=1.2,cex.axis=1.2,xlab="Urbanization", ylab= paste(input$attribute), type="p", pch=d$psmshape, col=dxy$cols)
    abline(v=Zcrit+.2, lty=2, lwd=2, col="blue")
    text(Zcrit,max(dxy$benefit),"zcrit", col= "blue")
    mtext(side=1,"high",line=4,adj=1,cex=0.8)
    mtext(side=1,"low",line=4,adj=0,cex=0.8)
    score_cohopres_m<-dxy
    legend("topleft", legend=c("Highest priority","Lowest priority"), pch=19,col=c(cols[1],cols[length(cols)]), cex=.8, bty="n")
  })
  
  output$mymap <- renderLeaflet({
    m<-leaflet(map_pred)%>%
      addTiles()%>%
      setView(lng=-120,lat=48,zoom=10)
  })
  
}

shinyApp(ui = ui, server = server)
