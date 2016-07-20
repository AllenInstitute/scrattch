library(shiny)
library(ggplot2)
library(dplyr)
library(Hmisc)
library(DBI)
library(RSQLite)
library(gridExtra)

shinyServer(function(input, output, session) {
  
  
  # URL Query Parsing
  
  ui_vals <- reactive({
    vals <- list(db = "//AIBSData2/mct-t200/ShinyApps/heatmapy/rpkm14.db",
                 genes = "Gad2,Slc17a6,Cux2",
                 primary = "cluster",
                 secondary = "cre",
                 clusters = "1:19",
                 scale = "scale.log",
                 showall = F,
                 sort = "none",
                 sortgene = "Gad2")
    
    if(length(session$clientData$url_search) > 0) {
      
      query <- as.list(parseQueryString(session$clientData$url_search))
      
      if(!is.null(query$db)) { vals$db = query$db }
      if(!is.null(query$genes)) { vals$genes = query$genes }
      if(!is.null(query$primary)) { vals$primary = query$primary }
      if(!is.null(query$secondary)) { vals$secondary = query$secondary }
      if(!is.null(query$clusters)) { vals$clusters = query$clusters }
      if(!is.null(query$scale)) { vals$scale = query$scale }
      if(!is.null(query$showall)) { vals$showall = as.logical(query$showall) }
      if(!is.null(query$sort)) { vals$sort = query$sort }
      if(!is.null(query$sortgene)) { vals$sortgene = query$sortgene }
      
    }
    
    return(vals)
  })
  
  # UI Elements with defaults set by URL query
  
  output$database_selection <- renderUI({
    textInput(inputId = "db", label = strong("Database Location"), value = ui_vals()$db, width = "100%")
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
  
  
  ################################################
  ## Primary and Secondary Annotation Selection ##
  ################################################
  
  desc <- reactive({
    con <- dbConnect(RSQLite::SQLite(),ui_vals()$db)
    get <- paste("SELECT * FROM desc",sep="")
    res <- dbSendQuery(con,get)
    desc <- dbFetch(res,n=-1)
    dbClearResult(res)
    
    return(desc)
  })
  
  annotations <- reactive({
    
    desc.list <- desc()$base
    names(desc.list) <- desc()$name
    
    return(desc.list)
  })

  output$primary_selection <- renderUI({
    selectInput("primary","Primary Grouping/Bottom Labels",annotations(),ui_vals()$primary)
  })
  
  output$secondary_selection <- renderUI({
    selectInput("secondary","Secondary Grouping/Top",annotations(),ui_vals()$secondary)
  })
  
  output$showall_selection <- renderUI({
    checkboxInput("showall","",ui_vals()$showall)
  })
  
  output$tertiary_selection <- renderUI({
    tert_annotations <- annotations()[!annotations() %in% c(input$primary,input$secondary)]
    selectizeInput("tertiary","Additional Annotations", choices = tert_annotations, multiple = T)
  })
  
  output$sortmode_selection <- renderUI({
    selectInput("sort",strong("Sort"),c("All"="all","Within Groups"="within","None"="none"),ui_vals()$sort)
  })
  
  output$scaling_selection <- renderUI({
    radioButtons("scale","Data Rescaling Options",
                 c("RPKM"="scale.abs",
                   "Log10(RPKM + 1)"="scale.log",
                   "RPKM/Max(RPKM) per Gene"="scale.rel",
                   "Log(RPKM+1)/Max(Log(RPKM+1))"="scale.log.rel"),
                 selected=ui_vals()$scale)
  })
  
  
  primary.anno <- reactive({
    
    # Select primary annotations
    primary <- list(id=paste0(input$primary,"_id"),
                    label=paste0(input$primary,"_label"),
                    color=paste0(input$primary,"_color")
                    )
    
    return(primary)
    
  })
  
  secondary.anno <- reactive({
    
    secondary <- list(id=paste0(input$secondary,"_id"),
                    label=paste0(input$secondary,"_label"),
                    color=paste0(input$secondary,"_color")
    )
    
    return(secondary)
    
  })
  
  pclusters <- reactive({
    clusters <- eval(parse(text=paste0("round(c(",input$clusters,"),0)")))
    
    primary <- primary.anno()
    con <- dbConnect(RSQLite::SQLite(),input$db)
    get <- paste0("SELECT ",primary$id," FROM anno")
    res <- dbSendQuery(con,get)
    all.clusters <- unique(dbFetch(res,n=-1))[,1]
    dbClearResult(res)
    
    intermediates <- all.clusters[round(floor(all.clusters),0) %in% clusters & all.clusters %% 1 != 0]
    # for the latest version of the plot, transitional cells have cluster ids with primary.secondary
    cluster_levels <- numeric()
    for(i in 1:length(clusters)) {
      cluster_levels <- c(cluster_levels,clusters[i],intermediates[round(floor(intermediates),0)==clusters[i]])
    }
    return(cluster_levels)
  })
  
  get_genes <- reactive({
    
    con <- dbConnect(RSQLite::SQLite(),input$db)
    get <- "SELECT gene FROM data"
    res <- dbSendQuery(con,get)
    all.genes <- dbFetch(res,n=-1)[,1]
    dbClearResult(res)
    dbDisconnect(con)
    
    genes <- unique(strsplit(input$genes,"[, \t]+")[[1]])
    
    genes <- genes[genes %in% all.genes]
    
    return(list(genes = genes))
    
  })
  
  get_data <- reactive({
    
    genes <- get_genes()$gene
    
    ##########################
    ## Annotation Retrieval ##
    ##########################
    
    con <- dbConnect(RSQLite::SQLite(),input$db)
    get <- paste("SELECT * FROM anno",sep="")
    res <- dbSendQuery(con,get)
    anno <- dbFetch(res,n=-1)
    dbClearResult(res)
    dbDisconnect(con)
    
    primary <- primary.anno()
    secondary <- secondary.anno()
    
    ## Check to see if all of the gene names are available in the database, and that there aren't blanks or other problems
    if(length(genes) > 0) {
    
      #######################
      ## Retrieval of data ##
      #######################
      
      ## Connect to the database and retrieve the rows for the selected genes.
      con <- dbConnect(RSQLite::SQLite(),input$db)
      getgenes <- paste("\"",genes,"\"",collapse=",",sep="")
      get <- paste("SELECT * from data WHERE gene IN (",getgenes,")",sep="")
      res <- dbSendQuery(con,get)
      genes.df <- data.frame(t(dbFetch(res)),stringsAsFactors=F)
      dbClearResult(res)
      dbDisconnect(con)
      
      ## Add an X to the beginning of gene names that start with a number
      ## This is necessary because the gene names will be used as column names, and R adds an X to colnames that start with numbers
      for(i in 1:length(genes)) {
        if (grepl("^[0-9]",genes[i])) {
          genes[i] <- paste("X",genes[i],sep="")
        }
        genes[i] <- sub("-",".",genes[i])
      }
      
      # Convert the gene names row to column names
      genes.cn <- sub("-",".",genes.df[1,])
      genes.df <- genes.df[2:nrow(genes.df),]
      
      # Convert the genes data frame to numeric
      genes.df <- sapply(genes.df,as.numeric)
      genes.df <- data.frame(genes.df)
      names(genes.df) <- genes.cn
      
      
      ###############
      ## Filtering ##
      ###############
      
      clusters <- pclusters()
      
      # Join the annotation and genes data frames
      sub.df <- data.frame(anno,genes.df)
      sub.df <- sub.df[sub.df[,primary$id] %in% clusters,]      
      
      #############
      ## Scaling ##
      #############
      
      ## If the y-axis is plotted on a log scale, add 1 to the data values to plot data + 1
      if(input$scale == "scale.log") {
        for(gene in genes) {
          sub.df[,gene] <- log10(sub.df[,gene] + 1)
        }
      }
      if(input$scale == "scale.rel") {
        for(gene in genes) {
          sub.df[,gene] <- sub.df[,gene]/max(sub.df[,gene])
        }
      }
      if(input$scale == "scale.log.rel") {
        for(gene in genes) {
          sub.df[,gene] <- log10(sub.df[,gene]+1)/log10((max(sub.df[,gene]+1)))
        }
      }
      
      #############
      ## Sorting ##
      #############
      
      if(input$sort == "all") {
        sort.df <- sub.df %>% arrange_(paste0("-",genes[1]),primary$id) %>% mutate(xpos=1:nrow(sub.df))
      } else {
        cluster_order <- data.frame(clust=clusters,
                                    plot_order=1:length(clusters))
        names(cluster_order)[1] <- primary$id
        sub.df <- sub.df %>% left_join(cluster_order,by=primary$id)
        
         # Sort the filtered data by cluster label followed by cre label, then add an x-position column
        if(input$sort == "within") {
          sort.df <- sub.df %>% 
            arrange_("plot_order",secondary$id,paste0("-",genes[1])) %>% 
            mutate(xpos=1:nrow(sub.df))  
        } else {
          # Sort the filtered data by primary label followed by secondary label, then add an x-position column
          sort.df <- sub.df %>% 
            arrange_("plot_order",secondary$id) %>% 
            mutate(xpos=1:nrow(sub.df))
        }
      }
      
      results <- list(data=sort.df,genes=genes)
      
      return(results)
      
    }
    
  })
  
  output$sort_selection <- renderUI({

    genes <- get_data()$genes
    genes <- sub("\\.","-",genes)
    selectInput("sortgene",strong("Sort By Gene"),genes,ui_vals()$sortgene)
    
  })
  
  ## Main ggpplot2 generating function
  buildplot <- function(pfontsize=14,expand=F) {
    data <- get_data()$data
    genes <- get_genes()$genes
    
    primary <- primary.anno()
    secondary <- secondary.anno()
    
    if(is.data.frame(data)) {
      
      sort.df <- data
      # For X-axis labels, we'll also need annotations for the mean of the positions of each cluster
      sort.lab <- sort.df %>% 
        group_by_(primary$id,primary$label) %>% 
        summarise(xmean=mean(c(min(xpos)-1,xpos)),y=length(genes)+1) %>%
        select_(primary$label,"xmean","y")
      
      colors <- colorRampPalette(c("darkblue","white","red"))(1001)
      
      if(input$autorange == "auto") {
        min.val <- 0
        max.val <- max(unlist(sort.df[,genes]))
      } else if (input$autorange == "manual") {
        min.val <- as.numeric(input$minrange)
        max.val <- as.numeric(input$maxrange)
      }
      
      ## Convert data to geom_rect() compatible table
      plot.df <- data.frame(xmin=numeric(),xmax=numeric(),ymin=numeric(),ymax=numeric(),fill=character())
      for(i in 1:length(genes)) {
        fill_ids <- round( (sort.df[,genes[i]] - min.val) / (max.val - min.val) * 1000 ) + 1
        fill_ids[fill_ids < 1] <- 1
        fill_ids[fill_ids > 1001] <- 1001
        gene.plot <- data.frame(xmin = sort.df$xpos - 1,
                                xmax = sort.df$xpos,
                                ymin = length(genes) - i,
                                ymax = length(genes) - i + 1,
                                fill = colors[fill_ids])
        plot.df <- rbind(plot.df,gene.plot)
      }
            
      primary.plot <- data.frame(xmin=sort.df$xpos-1,xmax=sort.df$xpos,ymin=-0.5,ymax=0,
                            fill=sort.df[,primary$color])
      secondary.plot <- data.frame(xmin=sort.df$xpos-1,xmax=sort.df$xpos,ymin=length(genes),ymax=length(genes)+0.5,
                            fill=sort.df[,secondary$color])
      plot.df <- rbind(plot.df,primary.plot,secondary.plot)
      
      ## add additional secondary color bars
      all.desc <- desc()
      primary.name <- all.desc$name[all.desc$base == input$primary]
      secondary.name <- all.desc$name[all.desc$base == input$secondary]
      other.desc <- all.desc[!all.desc$base %in% c(input$primary,input$secondary),]
      anno.color <- paste0(other.desc$base,"_color")
      anno_y_labels <- data.frame(breaks=numeric(),labels=character())
      
      if(input$showall) {
        for(j in 1:nrow(other.desc)) {
          anno.plot <- data.frame(xmin=sort.df$xpos-1,xmax=sort.df$xpos,
                                  ymin=length(genes)+j*0.5,ymax=length(genes)+j*0.5+0.5,
                                  fill=sort.df[,anno.color[j]])
          plot.df <- rbind(plot.df,anno.plot)
          
          anno_y <- data.frame(breaks=length(genes)+j*0.5+0.25,labels=other.desc$name[j])
          anno_y_labels <- rbind(anno_y_labels,anno_y)
        }
      } else {
      
        if(length(input$tertiary) > 0) {
          
          tert.desc <- other.desc[other.desc$base %in% input$tertiary,]
          
          for(j in 1:nrow(tert.desc)) {
            anno.plot <- data.frame(xmin=sort.df$xpos-1,xmax=sort.df$xpos,
                                    ymin=length(genes)+j*0.5,ymax=length(genes)+j*0.5+0.5,
                                    fill=sort.df[,paste0(input$tertiary[j],"_color")])
            plot.df <- rbind(plot.df,anno.plot)
            
            anno_y <- data.frame(breaks=length(genes)+j*0.5+0.25,labels=tert.desc$name[tert.desc$base == input$tertiary[j]])
            anno_y_labels <- rbind(anno_y_labels,anno_y)
          }
        }
      }
      
      ## build new, more complex y-axis labels
      y_labels <- data.frame(breaks=(1:length(genes)-0.5),labels=rev(genes))
      primary_y_label <- data.frame(breaks=-0.25,labels=primary.name)
      secondary_y_label <- data.frame(breaks=length(genes)+0.25,labels=secondary.name)
      y_labels <- rbind(y_labels,primary_y_label,secondary_y_label,anno_y_labels)
      
 #     hlines <- data.frame(yintercept=c(0:length(genes),length(genes)+0.5*(1:nrow(other.desc))))
      
      ##############
      ## Plotting ##
      ##############
      
      p <- ggplot() + geom_rect(data=plot.df,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill))
      
      p <- p + scale_fill_identity(guide=F)
      
  #    p <- p + geom_hline(data=hlines,aes(yintercept=yintercept))
      
      ##
      ## Axes and Cell labels
      ##
      
      # X-Axis
      p <- p + scale_x_continuous(breaks=sort.lab$xmean,labels=unlist(sort.lab[primary$label]),expand=c(0,0))
      # Y-Axis
      if(expand) {
       lab.len <- sort.df %>%
         select(one_of(input$cell.label)) %>% 
         rename_("lab"=input$cell.label) %>% 
         unique() %>%
         mutate(nchar=nchar(as.character(lab))) %>% 
         select(nchar) %>% max() %>% as.numeric()
       p <- p + geom_text(data=sort.df,aes_string(x="xpos-0.5",y=length(genes)+0.05*length(genes),label=input$cell.label,angle=90,hjust=0))
       p <- p + scale_y_continuous(breaks=y_labels$breaks,labels=y_labels$labels,limits=c(-0.5,length(genes)+0.028*lab.len*length(genes)),expand=c(0,0))     
      } else {
       p <- p + scale_y_continuous(breaks=y_labels$breaks,labels=y_labels$labels,expand=c(0,0))
      }
      
      ##
      ## Theme and other options
      ##
      
      p <- p + theme_classic(base_size=pfontsize)
      p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
      p <- p + theme(legend.title = element_blank())
      p <- p + theme(axis.title.y = element_blank(),axis.title.x = element_blank())

      ## Option to remove the legend entirely
#      if(!input$legend) {
#        p <- p + theme(legend.position="none")
#      }
      
      ##
      ## Output the plot
      ##
      
      return(p)

      
    } else {
      unmatched.genes <- data
      ## If the gene isn't found in the table, print an error message:
      err.df <- data.frame(x=c(1,1),y=c(1,0.75),lab=c(paste0(paste0(unmatched.genes,collapse=",")," can't be found in our database."),
                                                      "Check the spelling, or search NCBI Gene for synonyms in the Also Known As field."))
      p <- ggplot(err.df,aes(x=x,y=y)) +
        geom_text(aes(label=lab)) + ylim(0,1.25) + theme_classic()
      
      return(p)
    }

  }
  
  topwidth <- reactive({
    
#     if(input$expandtop) {
#       
#       # Retrieve cell annotations (anno) and the list of all gene names in the database
#       con <- dbConnect(RSQLite::SQLite(),input$db)
#       get <- paste("SELECT * FROM anno",sep="")
#       res <- dbSendQuery(con,get)
#       anno <- dbFetch(res,n=-1)
#       dbClearResult(res)
#       dbDisconnect(con)
#       
#       cats <- 1:7
#       
#       if(input$filters == "all") { cats <- 1:7 }
#       if(input$filters == "neuron") { cats <- c(1,2) }
#       if(input$filters == "glia") { cats <- 3:6 }
#       if(input$filters == "inhib") { cats <- 1 }
#       if(input$filters == "excit") { cats <- 2 }
#             
#       ncells <- sum(anno$cat_id %in% cats)
#       
#       w <- 150 + ncells * 15
#     
#     } else {
      w <- "auto"
#    }
    
    return(w)
    
  })
  
  topheight <- reactive({
    h <- as.numeric(input$plotheight)
    return(h)
  })
  
  output$ptop <- renderPlot({
    
    try(buildplot(pfontsize=as.numeric(input$plotfont)),silent=T)

  },width=topwidth,height=topheight)
  
  click_data <- reactive({
    
    data <- get_data()$data
    genes <- get_data()$genes
    annos <- annotations()
    
    colors <- colorRampPalette(c("darkblue","white","red"))(1001)
    max.val <- max(unlist(data[,genes]))
    
    if(length(input$plot_click) > 0) {
      x <- floor(input$plot_click$x)
      y <- floor(input$plot_click$y)
      
      clicked.cell.row <- data %>% filter(xpos==x+1)
      clicked.cell.anno <- clicked.cell.row %>% select(-one_of(genes))
      
      if(y %in% (1:length(genes)-1)) {
        clicked.gene <- genes[length(genes)-y]
        clicked.expression <- round(clicked.cell.row[,clicked.gene],2)
        clicked.expression.color <- colors[floor(clicked.expression/max.val*1000)+1]
      } else {
        clicked.gene <- "None Selected"
        clicked.expression <- "N/A"
        clicked.expression.color <- "black"
      }

      click.data <- data.frame(label=c("Cell ID",clicked.gene),
                               fill=c("white",clicked.expression.color),
                               value=c(clicked.cell.row$sample_id,clicked.expression),
                               y=c(0,-1.1))
      
      for(i in 1:length(annos)) {
        anno.data <- data.frame(label=names(annos)[i],
                                fill=clicked.cell.row[,paste0(annos[i],"_color")],
                                #value=paste(clicked.cell.row[,paste0(annos[i],"_id")],clicked.cell.row[,paste0(annos[i],"_label")]),
                                value=paste0(clicked.cell.row[,paste0(annos[i],"_label")]),
                                y=-i*1.1-1.1)
        click.data <- rbind(click.data,anno.data)
      }
    } else {
      click.data <- data.frame(label=c("Cell ID","None Selected"),
                               fill="white",
                               value=c("None Selected","N/A"),
                               y=c(0,-1.1))
      for(i in 1:length(annos)) {
        anno.data <- data.frame(label=names(annos)[i],
                                fill="white",
                                value="N/A",
                                y=-i*1.1-1.1)
        click.data <- rbind(click.data,anno.data)
      }
    }
    
    return(click.data)
    
  })

  output$click_plot <- renderPlot({
    try(  
    ggplot(click_data()) + 
      geom_rect(aes(xmin=0,xmax=1,ymin=y,ymax=y-1,fill=fill),color="black") +
      geom_text(aes(x=-0.1,y=y-0.5,label=label),hjust=1,size=4.25) +
      geom_text(aes(x=1.1,y=y-0.5,label=value),hjust=0,size=4.25) +
      scale_x_continuous(limits=c(-3,4),expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_identity() +
      theme_classic(5) +
      theme(axis.title=element_blank(),
            axis.text=element_blank(),
            axis.ticks=element_blank(),
            axis.line=element_blank())
    ,silent=T)  
  })


  build_legend_plot <- function(pfontsize = 14) {
    
    data <- get_data()$data
    genes <- get_data()$genes
    
    colors <- colorRampPalette(c("darkblue","white","red"))(1001)
    
    if(input$autorange == "auto") {
      min.val <- 0
      max.val <- max(unlist(data[,genes]))
    } else if (input$autorange == "manual") {
      min.val <- as.numeric(input$minrange)
      max.val <- as.numeric(input$maxrange)
    }
    
    ## Build geom_rect() compatible table
    legend_data <- data.frame(xmin = 1:1001,
                              xmax = 1:1001+1,
                              ymin = 0,
                              ymax = 1,
                              fill = colors)
    
    if(input$scale == "scale.abs") {
      scale_name <- "RPKM"
    } else if(input$scale == "scale.log") {
      scale_name <- "log10(RPKM + 1)"
    } else if(input$scale == "scale.rel") {
      scale_name <- "RPKM/max(RPKM)"
    } else if(input$scale == "scale.log.rel") {
      scale_name <- "log10(RPKM + 1)/max(log10(RPKM + 1))"
    }
    
    segment_data <- data.frame()
    
    legend_plot <- ggplot(legend_data) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
      geom_segment(aes(x = min(xmin), xend = max(xmax), y = 0, yend = 0)) +
      scale_fill_identity() +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(scale_name, breaks=c(0,250,500,750,1000),
                         labels=round(seq(min.val, max.val, by = (max.val-min.val)/4),2)) +
      theme_classic(base_size = pfontsize) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.x = element_blank())
    
    return(legend_plot)
  }  

  output$legend_plot <- renderPlot({
    
    build_legend_plot(as.numeric(input$plotfont))
    
  })


  output$downloadPlot <- downloadHandler(    
    
    filename = function() { paste("heatmap", '.pdf', sep='') },
    content = function(file) {
      heatmap_plot <- buildplot(pfontsize=as.numeric(input$dlf)) + theme(line=element_line(size=0.4))
      legend_plot <- build_legend_plot(pfontsize=as.numeric(input$dlf)) + theme(line=element_line(size=0.4),plot.margin = unit(c(0.1,0.35,0.1,0.35),"npc"))
      
      plot_list <- list(heatmap_plot,legend_plot)
      out_h <- as.numeric(input$dlh)
      out_w <- as.numeric(input$dlw)
      
      device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      ggsave(file, plot = arrangeGrob(grobs = plot_list,
                                      heights = c(out_h/8*7,out_h/8)),
             device = device, width=as.numeric(input$dlw), height=as.numeric(input$dlh))
    }
  )
  
  
  
  output$url <- renderUI({
    
    url <- "http://ibs-bosiljkat-ux1:8787/heatmap_it/?"
    
    vals <- list(db = input$db,
                 genes = gsub("[\t ,]+",",",input$genes),
                 primary = input$primary,
                 secondary = input$secondary,
                 clusters = input$clusters,
                 scale = input$scale,
                 showall = as.character(input$showall),
                 sort = input$sort,
                 sortgene = input$sortgene)
    
    for(i in 1:length(vals)) {
      url <- paste0(url,names(vals)[i],"=",vals[[i]])
      if(i < length(vals)) {
        url <- paste0(url,"&")
      }
    }
    
    a("DL",href=url)
  })
  
    
})
