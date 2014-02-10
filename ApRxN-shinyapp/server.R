library(shiny)
library(ggplot2)
library(stringr)
library(plyr)

# load data
gxdata <- read.csv("Ap.dt.csv")

# Define server logic required to plot variables 
shinyServer(function(input, output) {
  
  # Caption text
  captionText <- reactive({
    paste("Reaction norm for ", input$searchtext, "related genes")
  })
  
  # Return the formula text for printing as a caption
  output$caption <- renderText({
    captionText()
  })
  
  # Subset data to transcripts that have a match to 'searchtext' in best.hit.to.nr or GO.Biological.Process
  datasub <- reactive({
    sd <- gxdata[union(grep(input$searchtext, gxdata$best.hit.to.nr), grep(input$searchtext, gxdata$GO.Biological.Process)), ]
  })
  
  # Extract GO annotation to display
  output$GOannotation <- renderPrint({
    foo <- datasub()
    unlist(str_split(foo[1, "GO.Biological.Process"], " ", n=2))[2]
  })
  
  # Generate a plot of the requested variable against temperature
  output$RxNplot <- renderPlot({
    gdata <- datasub()
    #cat("gdata", dim(gdata), "\n")
    g <- ggplot(gdata, aes(x=val, y=exp.scaled, group=Transcript)) + 
      geom_smooth(aes(colour = exp_type)) + 
      facet_grid(. ~ colony)
    print(g)
  })
  
})