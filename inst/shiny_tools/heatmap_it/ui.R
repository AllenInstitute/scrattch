library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"),
  # Application title
  fluidRow(column(3,h2("Heatmap Generator")),
           column(8,uiOutput("database_selection")),
           column(1,uiOutput("url"))
  ),
  fluidRow(
           column(3,uiOutput("primary_selection")),
           column(3,uiOutput("secondary_selection")),
           column(4,
                  uiOutput("tertiary_selection")
                  
           ),
           column(2,
                  strong("Show All"),
                  uiOutput("showall_selection"))
  ),

  fluidRow(
    column(1,uiOutput("sortmode_selection")
    ),
    column(2,
           uiOutput("sort_selection")
           ),
    column(7,
           uiOutput("gene_selection"),
           tags$style(type='text/css', "#genes { width: 700px; }")
    ),
    column(2,
           uiOutput("cluster_selection")
    )
  ),
  
  plotOutput("ptop",height="100%",
             click=clickOpts(id="plot_click",clip=TRUE)),
  fluidRow(
    column(3,
           plotOutput("click_plot",height=400)
           ),
    column(3,
           uiOutput("scaling_selection"),
           strong("Range Options"),
           fluidRow(
             column(4,radioButtons("autorange","",
                                   c("Auto" = "auto",
                                     "Manual" = "manual"))),
             column(4,textInput("minrange","Min","0")),
             column(4,textInput("maxrange","Max","5"))
           )
    ),
    column(3,
           plotOutput("legend_plot",height=85),
      strong("Display Options"),
#      checkboxInput(inputId="expandtop",label="Label Cells",value=F),
#      checkboxInput(inputId="expandbottom",label="Label Bottom Cells",value=F),
#      checkboxInput(inputId="legend",
#                    label="Show Legend",
#                    value=T),
      fluidRow(
        column(6,textInput("plotheight","Height (px)","500")),
        column(6,textInput("plotfont","Font Size","14"))
        )
      
    ),
    column(3,
#           selectInput("db.version",label="",
#                       choices=c("Anderson VMH"),"Anderson VMH"),
           strong("Download Options"),
           fluidRow(
             column(4,textInput("dlw","Width (in)",8)),
             column(4,textInput("dlh","Height (in)",4)),
             column(4,textInput("dlf","Font (pt)",6))),
      downloadButton('downloadPlot')
    )
    
  )
  
))
