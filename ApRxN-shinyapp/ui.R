library(shiny)

# Define UI for selecting GO term
shinyUI(pageWithSidebar(
    
  # Application title
  headerPanel("Aphaenogaster thermally-responsive genes"),
  
  # Sidebar to select which GO term to plot transcripts for
  sidebarPanel(
    textInput("searchtext", "Search term:", value = "response to stress"),
    submitButton(text = "Apply Changes"),
    br(),
    p("Search terms to try: heat shock protein, behavior, phagocytosis, immune system, GO:0006950")
  ),
  
  # Display caption and plot of requested GO term
  mainPanel(
    h3(textOutput("caption")),
    
    plotOutput("RxNplot"),
    
    verbatimTextOutput("GOannotation")
    )
))
