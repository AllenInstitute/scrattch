library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  fluidRow(
    column(6,
           uiOutput("database_selection")
           ),
    column(5,
           uiOutput("primary_selection")
           ),
    column(1,
           uiOutput("url")
           )
  ),
  fluidRow(
    column(3,
           uiOutput("plot_type_selection")
           ), #End plot type column
    column(3,
           uiOutput("gene_selection")
           ), # End gene input column
    column(1,strong("Log Scale"),br(),
           uiOutput("log_selection")
      ), # End logscale column
    column(2,
           strong("Square Labels"),br(),
           uiOutput("square_selection")
    ),
    column(3,
           uiOutput("cluster_selection")
           ) # End cluster selection column
    ), # End options row
  fluidRow(
    plotOutput(outputId = "outplot",
               width =  "100%",
               height = "800px")
  ),# End plot row
  fluidRow(
    p(paste0("scrattch version ",packageVersion("scrattch")))
  )
  
) # End fluidPage
) # End shinyUI
