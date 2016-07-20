library(shiny)
library(scrattch)
library(DBI)
library(RSQLite)

shinyServer(function(input, output, session) {
  
  # URL Query Parsing
  
  ui_vals <- reactive({
    vals <- list(db = "/data/mct-t200/ShinyApps/heatmapy/rpkm14.db",
                 ptype = "heatcell",
                 genes = "Gad2,Slc17a6,Cux2",
                 group = "final",
                 clusters = "1:49",
                 log = T,
                 square = F)
    
    if(length(session$clientData$url_search) > 0) {
      
      query <- as.list(parseQueryString(session$clientData$url_search))
      
      if(!is.null(query$db)) { vals$db = query$db }
      if(!is.null(query$ptype)) { vals$ptype = query$ptype }
      if(!is.null(query$genes)) { vals$genes = query$genes }
      if(!is.null(query$group)) { vals$group = query$group }
      if(!is.null(query$clusters)) { vals$clusters = query$clusters }
      if(!is.null(query$log)) { vals$log = as.logical(query$log) }
      if(!is.null(query$square)) { vals$square = as.logical(query$square) }
      
      
    }
        
    return(vals)
  })
  
  # UI Elements with defaults set by URL query
  
  output$database_selection <- renderUI({
    textInput(inputId = "db", label = strong("Database Location"), value=ui_vals()$db,width="100%")
  })
  
  output$plot_type_selection <- renderUI({
    selectInput(inputId = "ptype",
                label = strong("Plot Type"),
                choices = c("Cell Heatmap" = "heatcell",
                            "Cell Barplot" = "barcell",
                            "Cluster Heatmap (mean)" = "heater.mean",
                            "Cluster Heatmap (median)" = "heater.median",
                            "Cluster Heatmap (trimmed mean)" = "heater.tmean",
                            "Cluster Heatmap (percent > 0)" = "heater.percent",
                            "Cluster Boxplot" = "boxter",
                            "Cluster Violins" = "pottery"),
                selected = ui_vals()$ptype
    )
  })
  
  output$gene_selection <- renderUI({
    textInput(inputId = "genes",
              label = strong("Genes"),
              value = ui_vals()$genes)
  })
  
  output$cluster_selection <- renderUI({
    textInput(inputId = "clusters",
              label = strong("Clusters"),
              value = ui_vals()$clusters)
  })
  
  output$log_selection <- renderUI({
    checkboxInput(inputId = "log",
                  label = "",
                  value = ui_vals()$log)
  })
  
  output$square_selection <- renderUI({
      checkboxInput(inputId = "square",
                    label = "",
                    value = ui_vals()$square)
  })
  
  # Description retrieval for generation of the Grouping selection.
  
  desc <- reactive({
    if(input$db == "internal") {
      desc <- v1_desc
    } else {
      con <- dbConnect(RSQLite::SQLite(),input$db)
      get <- paste("SELECT * FROM desc",sep="")
      res <- dbSendQuery(con,get)
      desc <- dbFetch(res,n=-1)
      dbClearResult(res)
    }
    return(desc)
  })
  
  annotations <- reactive({
    
    desc.list <- desc()$base
    names(desc.list) <- desc()$name
    
    return(desc.list)
  })
  
  # Primary selection will either auto-get the grouping as the first option in the
  # desc table, or use the grouping provided by the URL query
  output$primary_selection <- renderUI({
    if(length(ui_vals()$group) == 0) {
      selectInput("group","Group cells by",annotations(),annotations()[1])
    } else {
      selectInput("group","Group cells by",annotations(),ui_vals()$group)
    }
  })
  
  
  output$outplot <- renderPlot({
 #   genes <- fix_mouse_genes(split_cst(input$genes))
    genes <- split_cst(input$genes)
    clusters <- chr_to_num(input$clusters)
    
    if(input$square) {
      labeltype <- "square"
    } else {
      labeltype <- "poly"
    }
    
    if(input$ptype == "heatcell") {
      heatcell_plot(genes = genes,
                    clusters = clusters,
                    grouping = input$group,
                    font = 12,
                    logscale = input$log,
                    data_source = input$db,
                    labeltype = labeltype)
    } else if(input$ptype == "barcell") {
      barcell_plot(genes = genes,
                   clusters = clusters,
                   grouping = input$group,
                   font = 12,
                   logscale = input$log,
                   data_source = input$db) 
    } else if(input$ptype == "heater.mean") {
      heater_plot(genes = genes,
                  clusters = clusters,
                  grouping = input$group,
                  font = 12,
                  logscale = input$log,
                  data_source = input$db,
                  calculation = "mean")
    } else if(input$ptype == "heater.tmean") {
      heater_plot(genes = genes,
                  clusters = clusters,
                  grouping = input$group,
                  font = 12,
                  logscale = input$log,
                  data_source = input$db,
                  calculation = "trimmed_mean")
    } else if(input$ptype == "heater.median") {
      heater_plot(genes = genes,
                  clusters = clusters,
                  grouping = input$group,
                  font = 12,
                  logscale = input$log,
                  data_source = input$db,
                  calculation = "median")
    } else if(input$ptype == "heater.percent") {
      heater_plot(genes = genes,
                  clusters = clusters,
                  grouping = input$group,
                  font = 12,
                  logscale = input$log,
                  data_source = input$db,
                  calculation = "percent")
    } else if(input$ptype == "boxter") {
      boxter_plot(genes = genes,
                  clusters = clusters,
                  grouping = input$group,
                  font = 12,
                  logscale = input$log,
                  data_source = input$db)
    } else if(input$ptype == "pottery") {
      pottery_plot(genes = genes,
                   clusters = clusters,
                   grouping = input$group,
                   font = 12,
                   logscale = input$log,
                   data_source = input$db)
    }
  
  })
  
  output$url <- renderUI({
    
    url <- "http://ibs-bosiljkat-ux1:8787/scrattch_it/?"
    
    vals <- list(db = input$db,
                 ptype = input$ptype,
                 genes = gsub("[\t ,]+",",",input$genes),
                 group = input$group,
                 clusters = input$clusters,
                 log = as.character(input$log),
                 square = as.character(input$square))
    
    for(i in 1:length(vals)) {
      url <- paste0(url,names(vals)[i],"=",vals[[i]])
      if(i < length(vals)) {
        url <- paste0(url,"&")
      }
    }
    
    a("DL",href=url)
  })
  

})
