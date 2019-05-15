## libraries
library(dplyr)
library(shiny)

## Read in the data/model estimates and join
psm_pre<-read.table("../analysis/results/PSM_predictions.txt",header=TRUE)
psm_pre3<-psm_pre[1:51,]
#spawn<-read.csv("../data/spawner_data.csv", header=TRUE)


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
      selectInput('benefit', 'Benefit', names(iris)),
      
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
    #create a variable of 0/1 for if PSM is less than threshold
    psm_pre3$psmcol<-"darkred" 
    psm_pre3$psmcol[psm_pre3$p.psm.mean<input$psm_thresh]<-"lightgreen"
    Zcrit<-min(psm_pre$Z.mean[psm_pre$p.psm.mean>input$psm_thresh])
    
    ##Calculate difference between Z_mean and Zcrit (=deltaZ, or the change in Z required to get PSM to 40%)
    psm_pre3$Zcrit<-Zcrit
    psm_pre3$deltaZ<-psm_pre3$Zcrit-psm_pre3$Z.mean
    
    plot(psm_pre3$p.psm.mean, psm_pre3$Z.mean, pch=19,col=psm_pre3$psmcol, cex.lab=1.2,cex.axis=1.2,cex=2, ylab="Urbanization effections (Z)", xlab= "Mean Pre-Spawn Mortality")
    
    abline(v=input$psm_thresh, lty=2, lwd=2)
    
    text(input$psm_thresh+.02,min(psm_pre3$Z.mean),label="PSM threshold", cex=1.2)
    
    abline(h=Zcrit, lty=2, lwd=2, col="blue")
    
    polygon(c(input$psm_thresh,1,1,input$psm_thresh),c(Zcrit,Zcrit,max(psm_pre3$Z.mean)+.5,max(psm_pre3$Z.mean)+.5),
            col=adjustcolor("salmon",alpha.f=0.5),
            border=NA)
    
    text(.02,Zcrit+.04,label="Zcrit", col="blue",cex=1.2)
    
    
  })
  
}

shinyApp(ui = ui, server = server)
