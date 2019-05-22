## libraries
library(dplyr)
library(shiny)

## Read in the data/model estimates and join
psm_pre<-read.table("../analysis/results/PSM_predictions.txt",header=TRUE)
spawn<-read.csv("../data/spawner_data.csv", header=TRUE)
spatial<-read.csv("../data/spatial_data.csv", header=TRUE)
#flags
allsites=FALSE#if choose FALSE, then will only use the 51 sites for which there are observations of PSM


#choose what data you want to include
allsites=FALSE #if false, selects out only sites with PSM calculated from field data, rather than sites with predicted PSM too

# 
ui <- fluidPage(
  
  # App title ----
  titlePanel("Prioritizing Coho Sites"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "psm_thresh",
                  label = "Pre-Spawn Mortality (PSM) Threshold",
                  min = 0,
                  max = 1.0,
                  value = 0.25)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Plot PSM vs. Z ----
      plotOutput(outputId = "plot_psm_z")
      
    )
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
    
    
    source ("../analysis/source/prepforplots.R")
    
    ##Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    d$Zcrit<-Zcrit
    d$deltaZ<-d$Zcrit-d$Z.mean
    
    plot(d$p.psm.mean, d$Z.mean, pch=19,col=d$psmcol, cex.lab=1.2,cex.axis=1.2,cex=2, ylab="Urbanization effects (Z)", xlab= "Mean Pre-Spawn Mortality")
    
    abline(v=input$psm_thresh, lty=2, lwd=2)
    
    text(input$psm_thresh+.02,min(d$Z.mean),label="PSM threshold", cex=1.2)
    
    abline(h=Zcrit, lty=2, lwd=2, col="blue")
    
    polygon(c(input$psm_thresh,1,1,input$psm_thresh),c(Zcrit,Zcrit,max(d$Z.mean)+.5,max(d$Z.mean)+.5),
            col=adjustcolor("salmon",alpha.f=0.5),
            border=NA)
    
    text(.02,Zcrit+.04,label="Zcrit", col="blue",cex=1.2)
    
    
  })
  
}

shinyApp(ui = ui, server = server)
