
## libraries
library(here)
library(dplyr)
library(RColorBrewer)
library(colorRamps)
## Step 1:  Read in the data/model estimates- use only predicted attributes for now
psm_pre <- read.table(here("analysis","results","PSM_predictions.txt"), header=TRUE)
spawn <- read.csv(here("data","spawner_data.csv"), header=TRUE)
spatial_pred <- read.csv(here("data","spatial_data_predict.csv"), header=TRUE)
#spatial<- read.csv(here("data","spatial_data.csv"), header=TRUE)

salmon <- read.csv(here("data","salmonscape","WA_integrated_Fish_coho_chinook_chum_with_WADOE.csv"),header=TRUE, skip=1)

## Step 2: Choices: select the threshold psm and you want to use, and select all IDs or only IDs for which we have PSM data (rather than predicted PSM)
input <- as.data.frame(NA)
input$psm_thresh <- 0.25
input$attribute<-"Coho_Presence_m"

predsites <- TRUE #if false, selects out only IDs with PSM calculated from field data, rather than IDs with predicted PSM

## Step 3: combine all the data and prep for plotting calculate mean spawner abundance by ID, across years
## combined data file with things we want to plot is called "d"
source(here("analysis","source","prepforplots.R"))
focalcols<-c("Coho_Presence_m","nsp_pres","ChinFa_Presence_m")

ui <- pageWithSidebar(
  
  # App title ----
 headerPanel("Prioritizing Coho Sites"),
    # Sidebar panel for inputs ----
 sidebarPanel(
      # Input: Selector for attribute to plot against z ----
      selectInput(inputId = "attribute", 
                  label = "Attribute:", 
                  focalcols,
                  selected=focalcols[1]),
      
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
      
      # Output: Plot selected attribute vs Z
      plotOutput(outputId = "plot_attribute")
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
    d$psmshape<-17 
    d$psmshape[d$p_psm_mean<input$psm_thresh]<-19
    ## Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    ## for all sites and select out just the bad sites
    #plot(dsubs$p_psm_mean, dsubs$Z_mean, pch=d$shape,cex.lab=1.2,cex.axis=1.2,cex=2, ylab="Urbanization effects (Z)", xlab= "Mean Pre-Spawn Mortality")
    plot(d$Z_mean,d$p_psm_mean,pch=d$psmshape, col="gray", cex.lab=1.2,cex.axis=1.2,cex=1.52, 
         xlab="Urbanization score (Z)", ylab= "Pre-Spawn Mortality")
    abline(h=input$psm_thresh, lty=2, lwd=2)
    text(min(d$Z_mean, na.rm=TRUE)+1,input$psm_thresh+.02,label="PSM threshold", cex=1.2)
    abline(h=Zcrit, lty=2, lwd=2, col="blue")
    
    polygon(c(input$psm_thresh,1,1,input$psm_thresh),c(Zcrit,Zcrit,max(d$Z.mean)+.5,max(d$Z.mean)+.5),
            col=adjustcolor("salmon",alpha.f=0.5),
            border=NA)
    
    text(.02,Zcrit+.04,label="Zcrit", col="blue",cex=1.2)
    
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
  
}

shinyApp(ui = ui, server = server)
