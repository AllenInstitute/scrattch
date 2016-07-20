library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  fluidRow(
    column(5,
           uiOutput("database_selection")
           ),
    column(3,
           uiOutput("group_selection")
           ),
    column(3,
           uiOutput("data_selection")
    ),
    column(1,
           uiOutput("url")
           )
  ),
  fluidRow(
    column(7,strong("Log Scale"),br(),
           uiOutput("log_selection")
      ), # End logscale column
    column(1,
           uiOutput("div_logic")
    ),
    column(2,
           uiOutput("div_selection")
    ),
    column(2,
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
