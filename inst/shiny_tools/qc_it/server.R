library(shiny)
library(dplyr)
library(scrattch)
library(ggplot2)
library(DBI)
library(RSQLite)
library(colorspace)

source("stats_functions.R")

shinyServer(function(input, output, session) {
  
  # URL Query Parsing
  
  ui_vals <- reactive({
    # set defaults
    vals <- list(db = "/data/rnaseqanalysis/SMARTer/SC_Core/analysis_1602/mouse_lgn_20160419.db",
                 group = "",
                 metric = "",
                 clusters = "1:49",
                 log = F,
                 div_logic = F,
                 div = "")
    
    if(length(session$clientData$url_search) > 0) {
      
      query <- as.list(parseQueryString(session$clientData$url_search))
      
      for(val in names(vals)) {
        if(!is.null(query[[val]])) {
          vals[[val]] <- query[[val]]
        }
      }
      
    }
    
    return(vals)
  })
  
  # UI Elements with defaults set by URL query
  
  output$database_selection <- renderUI({
    textInput(inputId = "db", label = strong("Database Location"), value=ui_vals()$db,width="100%")
  })
  
  output$cluster_selection <- renderUI({
    textInput(inputId = "clusters",
              label = strong("Clusters"),
              value = ui_vals()$clusters)
  })
  
  output$log_selection <- renderUI({
    checkboxInput(inputId = "log",
                  label = "",
                  value = as.logical(ui_vals()$log))
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
  
  anno <- reactive({
    if(input$db == "internal") {
      desc <- v1_anno
    } else {
      con <- dbConnect(RSQLite::SQLite(),input$db)
      get <- paste("SELECT * FROM anno",sep="")
      res <- dbSendQuery(con,get)
      anno <- dbFetch(res,n=-1)
      dbClearResult(res)
    }
    return(anno)
  })
  
  cat_annotations <- reactive({
    
    desc_all <- desc()
    
    desc.list <- list()
    
    for(i in 1:nrow(desc_all)) {
      
      desc_label <- paste0(desc_all$base[i],"_label")
      
      if(try(sum(is.na(as.numeric(anno()[,desc_label]))) > 0,silent = T)) {
        desc.list <- c(desc.list,desc_all$base[i])
        names(desc.list)[length(desc.list)] <- desc_all$name[i]
      }
      
    }
    
    return(desc.list)
  })
  
  num_annotations <- reactive({
    
    desc_all <- desc()
    
    desc.list <- character()
    
    for(i in 1:nrow(desc_all)) {
      
      desc_label <- paste0(desc_all$base[i],"_label")
      
      if(try(sum(is.na(as.numeric(anno()[,desc_label]))) == 0,silent = T)) {
        desc.list <- c(desc.list,desc_all$base[i])
        names(desc.list)[length(desc.list)] <- desc_all$name[i]
      }
      
    }
    
    return(desc.list)
  })
  
  # Primary selection will either auto-get the grouping as the first option in the
  # desc table, or use the grouping provided by the URL query
  output$group_selection <- renderUI({
    if(length(ui_vals()$group) == 0) {
      selectInput("group","Group cells by",cat_annotations(),cat_annotations()[1])
    } else {
      selectInput("group","Group cells by",cat_annotations(),ui_vals()$group)
    }
  })
  
  output$data_selection <- renderUI({
    if(length(ui_vals()$group) == 0) {
      selectInput("data","QC Metric",num_annotations(),num_annotations()[1])
    } else {
      selectInput("data","QC Metric",num_annotations(),ui_vals()$data)
    }
  })
  
  output$div_logic <- renderUI({
    checkboxInput("div_logic","Use Division?",as.logical(ui_vals()$div_logic))
  })
  
  output$div_selection <- renderUI({
    if(length(ui_vals()$group) == 0) {
      selectInput("div","Divide By",num_annotations(),num_annotations()[2])
    } else {
      selectInput("div","Divide By",num_annotations(),ui_vals()$div)
    }
  })
  
  
  output$outplot <- renderPlot({

    clusters <- chr_to_num(input$clusters)
    
    group_id <- paste0(input$group,"_id")
    group_label <- paste0(input$group,"_label")
    group_color <- paste0(input$group,"_color")
    
    data_col <- paste0(input$data,"_label")
    div_col <- paste0(input$div,"_label")
    
    group_order <- data.frame(clusters = clusters) %>%
      mutate(order = 1:n())
    names(group_order)[1] <- group_id
    
    filter_this <- paste0(group_id," %in% c(",input$clusters,")")

    anno_groups <- anno() %>%
      filter_(filter_this) %>%
      select(one_of(group_id,group_label,group_color)) %>%
      unique() %>%
      left_join(group_order)  %>%
      mutate(textcolor = "black")
    
    anno_groups$textcolor[as.data.frame(coords(as(hex2RGB(unlist(anno_groups[,group_color])),"polarLUV")))$L < 50] <- "white"
    
    plot_anno <- anno() %>%
      filter_(filter_this) %>%
      left_join(group_order)
    
    plot_anno[,data_col] <- as.numeric(plot_anno[,data_col])
    plot_anno[,div_col] <- as.numeric(plot_anno[,div_col])
    
    if(input$log) {
      plot_anno[,data_col] <- log10(plot_anno[,data_col])
      plot_anno[,div_col] <- log10(plot_anno[,div_col])
    }
    
    if(input$div_logic) {
      div_mutate <- paste0(data_col,"/",div_col)
      new_data <- plot_anno[,data_col]/plot_anno[,div_col]
      plot_anno[,data_col] <- new_data
      
    }
    
    stats <- summarySE(plot_anno, measurevar = data_col, groupvars = group_label) %>%
      mutate(median_lab = ifelse(median < 1,round(median,3),median)) %>%
      mutate(median_lab = ifelse(median >= 1 & median < 10, round(median,2),median_lab)) %>%
      mutate(median_lab = ifelse(median > 10,round(median,1),median_lab)) %>% 
      left_join(anno_groups)
    
    
    y.axis.line <- data.frame(x=0.5,xend=0.5,y=0,yend=max(stats$max))
    
    med_chars <- max(nchar(as.character(stats$median_lab)))
    label_chars <- max(apply(stats[,group_label],1,nchar))
    
    med_space <- med_chars/sum(med_chars,label_chars)
    label_space <- label_chars/sum(med_chars,label_chars)
    
    both_label_height <- 0.4*max(stats$max)
    
    med_height <- both_label_height*med_space
    
    median.panel <- data.frame(xmin=0.5,xmax=nrow(stats)+0.5,ymin=-med_height,ymax=0)
    
    p <- ggplot() +
      geom_jitter(data = plot_anno,aes_string(x="order",y=data_col),color="skyblue",position = position_jitter(width = .3,height = 0),size=0.8) +
      geom_errorbar(data=stats,aes(x=order,ymin=q25,ymax=q75),width=0.5,size=0.2) +
      geom_point(data=stats,aes(x=order,y=median),size=1.3,color="red") +
      # Cluster color panels
      geom_rect(data=stats,aes_string(xmin="order-0.5",xmax="order+0.5",ymin=-both_label_height,ymax=0,fill=group_color)) +
      # White background for medians
      geom_rect(data=median.panel,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="white"),alpha = 0.5) +
      # Median label text
      geom_text(data=stats,aes(x=order,y=-med_height*0.1,label=median_lab),angle=90,vjust=0.35,hjust=1) +
      # Cluster/Group text
      geom_text(data=stats,aes_string(x="order",y="-med_height - med_height*0.2",label=group_label,color="textcolor"),angle=90,hjust=1,vjust=0.35) +
      # Value = 0 baseline
      geom_hline(aes(yintercept=0),size=0.2) +
      # Median/Group separator line
      geom_hline(aes(yintercept=-med_height),size=0.2) +
      # Y axis line
      geom_segment(data=y.axis.line,aes(x=x,xend=xend,y=y,yend=yend),size=0.2) +
      scale_x_continuous("",expand=c(0,0),breaks=1:(length(clusters))) +
      scale_y_continuous(desc()$name[desc()$base == input$data],
                         expand=c(0,0),
                         limits=c(-both_label_height,max(stats$max)*1.05)) +
      scale_fill_identity() +
      scale_color_identity() +
      theme_bw(14) +
      theme(axis.text.x = element_blank()) +
      theme(axis.ticks.x = element_blank()) +
      theme(panel.border = element_blank()) +
      theme(line=element_line(size=0.2))
      
    return(p)
  })
  
  output$url <- renderUI({
    
    url <- "http://ibs-bosiljkat-ux1:8787/qc_it/?"
    
    vals <- list(db = input$db,
                 group = input$group,
                 metric = input$metric,
                 clusters = input$clusters,
                 log = as.character(input$log),
                 div_logic = as.character(input$log),
                 div = input$div)

    for(i in 1:length(vals)) {
      url <- paste0(url,names(vals)[i],"=",vals[[i]])
      if(i < length(vals)) {
        url <- paste0(url,"&")
      }
    }
    
    a("DL",href=url)
  })
  

})
