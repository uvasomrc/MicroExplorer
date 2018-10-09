#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# load libraries
library(shiny)
library(shinythemes)
library(tidyverse)
library(plotly)
library(biomformat)
library(magrittr)
library(reshape2)
library(vegan)

# import utility functions
source("./utils.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
   # Shiny theme
   theme = shinytheme("united"),
   
   # Application title
   titlePanel("MicroExplorer"),
   
   # Navbar Pages
   navbarPage("MicroExplorer",
    tabPanel("Upload Data",
      sidebarLayout(
        sidebarPanel(
          tags$div(tags$h3("Upload Files"),
                   style="margin-bottom:50px"),
          radioButtons("fileFormat", label="Select File Format",
                       choices=c("CSV","BIOM"), inline = TRUE),
          conditionalPanel(
            condition = "input.fileFormat == 'CSV'",
            fileInput("countFile", label="Choose Count Data CSV"),
            fileInput("taxaFile", label="Choose Taxonomy Data CSV"),
            fileInput("sampleFile", label="Choose Sample Data CSV")),
          conditionalPanel(
            condition = "input.fileFormat == 'BIOM'",
            fileInput("biomFile", label="Choose BIOM File")),
          actionButton("validate", "Load Data")),  
        mainPanel(
          textOutput("validMsg")
        )
      )
    ),
    tabPanel("Filter Data",
      sidebarLayout(
        sidebarPanel(
          uiOutput("filterUI")
        ),
        mainPanel(
          plotOutput("hist_valid"),
          plotOutput("hist_filtered")
        )
      )       
    ),
    tabPanel("Stacked Bars",
      sidebarLayout(
        sidebarPanel(
          uiOutput("sbUI")
        ),
        mainPanel(
          plotlyOutput("stackedbar")
        )
      )
    )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #############
  # Upload Data
  #############
  # Reactive Raw Data 
  rawData <- reactive({
    if (input$fileFormat == "CSV") {
      req(!is.null(input$countFile) & !is.null(input$taxaFile) & !is.null(input$sampleFile))
      csv2dat(input$countFile$datapath, input$taxaFile$datapath, input$sampleFile$datapath)
    } else if (input$fileFormat == "BIOM") {
      req(!is.null(input$biomeFile))
      biom2dat(input$biomFile$datapath)
    }
  })
   
  # Reactive Validate Data 
  validData <- eventReactive(input$validate, {
    validateInputs(rawData()$countData, rawData()$taxaData, rawData()$sampleData)  
  })
  # valid message
  output$validMsg <- renderText({
    req(!is.null(validData()))
    validData()$msg
  })
  
  
  #############
  # Filter Data
  #############
  # render Filtering UI
  output$filterUI <- renderUI({
    req(!(grepl("ERROR", validData()$msg)))
    list(
      tags$div(tags$h3("Filter Data"), style="margin-bottom:50px"),
      tags$div(tags$h4("Filter Samples by Sequencing Depth: ")),
      numericInput("seqDepth", "Min reads", 
                   min = 0, max = max(colSums(validData()$countData)),
                   value = 5000),
      tags$div(tags$h4("Filter Taxa by Prevalence: "), style="margin-top:50px"),
      numericInput("minAbund", "Minimum Abundance (%)", 
                   min = 0, max = 100, value = 0.01),
      numericInput("minSamples", "In at least N (%) Samples", 
                   min = 0, max = 100, value = 5),
      tags$div(tags$h4("Filter Samples by Metadata: "), style="margin-top:50px"),
      actionButton("filter", "Filter")
    )
  })
  
  # Filter Data
  filteredData <- eventReactive(input$filter, {
    filterData(validData()$countData, validData()$taxaData, validData()$sampleData,
               input$seqDepth, input$minAbund, input$minSamples)  
  })
  
  # filtered data main panel
  output$hist_valid <- renderPlot({
   req(!is.null(validData()) & !(grepl("ERROR", validData()$msg)))
   seqDepth <- colSums(validData()$countData)
   hist(seqDepth, breaks=20)
  })
  output$hist_filtered <- renderPlot({
    req(!is.null(filteredData()))
    seqDepth <- colSums(filteredData()$countData)
    hist(seqDepth, breaks=20)
  })
  
  
  ###################
  # Stacked Bar Plots
  ###################
  # render stacked bar UI
  output$sbUI <- renderUI({
    req(!is.null(filteredData()$countData))
    list(
      radioButtons("plotDataset", label=h4("Dataset"),
                   choices=c("Raw","Filtered"), inline = TRUE,
                   selected = "Filtered"),
      selectInput("taxaLevel", h4("Taxonomy Level"),
                  choices=base::colnames(filteredData()$taxaData)),
      radioButtons("taxa2Plot", label = h4("Taxas to plot"),
                   choices=c("All", "MostAbundant"), inline = TRUE),
      conditionalPanel(
        condition = "input.taxa2Plot == 'MostAbundant'",
        sliderInput("numTaxa2Plot", label = h5("Maximum Taxa to plot:"),
                    min = 1, max = 20, value = 10)
      ),
      radioButtons("sortMethod", label= h4("Sort Samples"),
                   choices = c("Decreasing Taxa Abunance", "Cluster by Dissimilarity")),
      selectInput("facetField", h4("Facet By"),
                  choices=c("None", base::colnames(filteredData()$sampleData)), selected = "None"),
      actionButton("sbplot", "Plot")
    )
  })
  
  # render stacked bar plot
  stackedBarPlot <- eventReactive(input$sbplot, {
    if (input$plotDataset == "Filtered") {
      myStackedBarPlot(filteredData()$countData, filteredData()$taxaData, filteredData()$sampleData,
                       input$taxaLevel, input$taxa2Plot, input$numTaxa2Plot, 
                       input$sortMethod, input$facetField) 
    }
  })
  output$stackedbar <- renderPlotly({
    req(!is.null(stackedBarPlot()))
    stackedBarPlot()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

