#' Barplots of gene expression of individual cells
#' 
#' This function will generate plots similar to those in Figure 3a-c of Tasic, et al. (2015).
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' See Examples for sensible PDF output options.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param data_source A character object defining where the data is stored. Currently only works with "internal"
#' @param sort Logical object, determines if cells will be sorted within their clusters (this means that each column will no longer represent a single cell)
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#' 
#' @return a ggplot2 plot object
#' 
#' @examples
#' barcell_plot()
#' 
#' my_genes <- c("Ercc6","Ercc8","Trp53","Pgbd5")
#' my_clusters <- c(1,5,9,10,24,37)
#' my_barcell_plot <- barcell_plot(my_genes,my_clusters,sort=T)
#' 
#' ggsave("plot_output.pdf",my_barcell_plot,height=0.2*length(my_genes)+2,width=4)
barcell_plot <- function(genes=c("Hspa8","Snap25","Gad2","Slc17a6"),grouping="final",clusters=1:49,
                         data_source="internal",sort=F,logscale=F,
                         fontsize=7,labelheight=25) {
  
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    data <- scrattch::v1_data
    all.anno <- scrattch::v1_anno
    
    data <- data %>%
      filter(gene %in% genes)
    
    # Reformat the retrieved data
    row.names(data) <- data[,1]
    data <- data %>% 
      select(-1) %>% 
      t() %>% 
      as.data.frame()
    
    data <- data %>%
      mutate(sample_id=row.names(data)) %>%
      select(one_of(c("sample_id",genes)))
    
    # Fix Ndnf
    genes[genes == "A930038C07Rik"] <- "Ndnf"
    names(data)[names(data) == "A930038C07Rik"] <- "Ndnf"
    
    genes[genes == "9630013A20Rik"] <- "X9630013A20Rik"
    names(data)[names(data) == "9630013A20Rik"] <- "X9630013A20Rik"
    
  } else {
    
    library(DBI)
    library(RSQLite)
    
    con <- dbConnect(RSQLite::SQLite(),data_source)
    get <- "SELECT * FROM anno;"
    res <- dbSendQuery(con,get)
    all.anno <- dbFetch(res,n=-1)
    dbClearResult(res)
    getgenes <- paste("\"",genes,"\"",collapse=",",sep="")
    get <- paste("SELECT * from data WHERE gene IN (",getgenes,")",sep="")
    res <- dbSendQuery(con,get)
    data <- dbFetch(res,n=-1)
    row.names(data) <- data[,1]
    data <- data %>% 
      select(-1) %>% 
      t() %>% 
      as.data.frame()
    
    data <- data %>%
      mutate(sample_id=row.names(data)) %>%
      select(one_of(c("sample_id",genes)))
  }
    
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    rename_("plot_id" = paste0(grouping,"_id"),
            "plot_label" = paste0(grouping,"_label"),
            "plot_color" = paste0(grouping,"_color")) %>%
    filter(plot_id %in% clusters) %>%
    arrange(plot_id) %>%
    mutate(xpos=1:n())
  
  # Calculate the number of clusters
  nclust <- length(unique(data$plot_id))
  
  #Calculate the height of the label:
  labheight <- length(genes)*(labelheight/100)/(1-labelheight/100)
  
  # Build cell type label polygons
  poly.data <- data %>% 
    group_by(plot_id) %>%
    summarise(x2=min(xpos)-1,
              x1=max(xpos),
              color=plot_color[1])
  poly.data <- poly.data %>%
    mutate(x3=nrow(data)*(1:nclust-1)/nclust,
           x4=nrow(data)*(1:nclust)/nclust,
           y4=length(genes)+1+labheight*0.1,y3=length(genes)+1+labheight*0.1,y2=length(genes)+1,y1=length(genes)+1)
  poly <- data.frame(id=rep(poly.data$plot_id,each=4),color=rep(poly.data$color,each=4))
  poly.x <- numeric()
  poly.y <- numeric()
  for(i in 1:nrow(poly.data)) {
    poly.x <- c(poly.x,poly.data$x1[i],poly.data$x2[i],poly.data$x3[i],poly.data$x4[i])
    poly.y <- c(poly.y,poly.data$y1[i],poly.data$y2[i],poly.data$y3[i],poly.data$y4[i])
  }
  poly <- cbind(poly,poly.x=poly.x,poly.y=poly.y)
  
  # Build the cell type label rectangles
  xlab.rect <- data.frame(xmin=poly.data$x3,
                          xmax=poly.data$x4,
                          ymin=poly.data$y3,
                          ymax=poly.data$y3+labheight*0.9,
                          color=poly.data$color,
                          label=unique(data$plot_label))
  
  # Calculate Plot Scale bars
  scale.bars <- data.frame(gene=genes) %>% 
    mutate(ymin=1:n()) %>%
    mutate(ymax=ymin+0.9,ymid=ymin+0.45) %>% mutate(xmin=-nrow(data)*0.01,xmax=-1)
  
  # Build the maximum value labels for the right edge
  max.vals <- data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% unlist()
  max.labels <- data.frame(x=nrow(data)+0.01*nrow(data),y=1:length(genes)+0.5,
                           label=sci_label(max.vals))
  max.header <- data.frame(x=nrow(data)+0.075*nrow(data)/2,y=length(genes)+1.5,label="Max data")
  
  # Scale the datas
  data <- data
  for(i in 1:length(genes)) {
    gene <- genes[i]
    if(logscale) {
      data[gene] <- log10(data[gene]+1)
      data[gene] <- data[gene]/max(data[gene])*0.9 + i
    } else {
      data[gene] <- data[gene]/max(data[gene])*0.9 + i
    }
  }
  
  background_data <- data.frame(xmin=0,xmax=max(data$xpos),ymin=1,ymax=length(genes)+1,fill="#ADCFE0")
  
  # Plot setup
  p <- ggplot(data) +
    scale_fill_identity() +
    theme_classic(base_size=fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    geom_rect(data=background_data,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill))
  

  p <- p + geom_segment(data=poly.data,aes(x=x1,xend=x1,y=1,yend=y1),size=0.2,color="gray60",linetype="dashed")
  
  
  # plot the bars for each gene
  for(i in 1:length(genes)) {
    # if sort is true, arrange the values within each row from high to low.
    if(sort) {
      data <- data %>% arrange_("plot_id",paste0("-",genes[i])) %>% mutate(xpos=1:n())
    }
    
    # plot the rectangles for the barplots
    p <- p + geom_rect(data=data,aes_string(xmin="xpos-1",xmax="xpos",ymin=i,ymax=genes[i],fill="plot_color"))
    
  }
  
  p <- p + 
    # Cluster labels at the top of the plot
    geom_rect(data=xlab.rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=color)) +
    geom_text(data=xlab.rect,aes(x=(xmin+xmax)/2,y=ymin+0.05,label=label),angle=90,vjust=0.35,hjust=0,size=pt2mm(fontsize)) +
    geom_polygon(data=poly,aes(x=poly.x,y=poly.y,fill=color,group=id)) +
    # Scale bar elements
    geom_hline(data=scale.bars,aes(yintercept=ymin),size=0.2) +
    geom_segment(data=scale.bars,aes(x=xmin,xend=xmax,y=ymid,yend=ymid),size=0.2) +
    geom_segment(data=scale.bars,aes(x=xmin,xend=xmax,y=ymax,yend=ymax),size=0.2) +
    geom_segment(data=scale.bars,aes(x=xmax,xend=xmax,y=ymin,yend=ymax),size=0.2) +
    # Maximum value labels at the right edge of the plot
    geom_text(data=max.header,aes(x=x,y=y,label=label),angle=90,hjust=0,vjust=0.5,size=pt2mm(fontsize)) +
    geom_text(data=max.labels,aes(x=x,y=y,label=label),hjust=0,vjust=0.5,size=pt2mm(fontsize))
  
  # Axis scales
  p <- p + scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),breaks=1:length(genes)+0.5,labels=genes)
  
  return(p)
  
}

#' Heatmaps of gene expression of individual cells
#' 
#' This function will generate plots similar to those in the shiny heatmap geneterator.
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param data_source A character object defining where the data is stored. Currently only works with "internal"
#' @param normalize_rows Logical object, determines if data are normalized to the maximum value for each gene. If FALSE, the heatmap is normalized to the maximum value across all genes.
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#' 
#' @return a ggplot2 plot object
#' 
#' @examples
#' heatcell_plot()
#' 
#' my_genes <- c("Ercc6","Ercc8","Trp53","Pgbd5")
#' my_clusters <- c(1,5,9,10,24,37)
#' heatcell_plot <- manybar_plot(my_genes,my_clusters,norm=T,font=12)
#' 
#' ggsave("plot_output.pdf",my_manybar_plot,height=0.2*length(my_genes)+2,width=4)
heatcell_plot <- function(genes=c("Hspa8","Snap25","Gad2","Slc17a6"),clusters=1:49,
                          data_source="internal",
                          sort=F,logscale=T,normalize_rows=F,
                          fontsize=7,labelheight=25) {
  
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    data <- scrattch::v1_data
    all.anno <- scrattch::v1_anno
    
    data <- data %>%
      filter(gene %in% genes)
    
    # Reformat the retrieved data
    row.names(data) <- data[,1]
    data <- data %>% 
      select(-1) %>% 
      t() %>% 
      as.data.frame()
    
    data <- data %>%
      mutate(sample_id=row.names(data)) %>%
      select(one_of(c("sample_id",genes)))
    
  }
  
  # Fix Ndnf
  genes[genes == "A930038C07Rik"] <- "Ndnf"
  names(data)[names(data) == "A930038C07Rik"] <- "Ndnf"
  
  genes[genes == "9630013A20Rik"] <- "X9630013A20Rik"
  names(data)[names(data) == "9630013A20Rik"] <- "X9630013A20Rik"
  
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    filter(final_id %in% clusters) %>%
    arrange(final_id) %>%
    mutate(xpos=1:n())
  
  #Calculate the height of the label:
  labheight <- length(genes)*(labelheight/100)/(1-labelheight/100)
  
  # Build cell type label polygons
  poly.data <- data %>% 
    group_by(final_id) %>%
    summarise(x2=min(xpos)-1,x1=max(xpos),color=final_color[1])
  poly.data <- poly.data %>%
    mutate(x3=nrow(data)*(1:length(clusters)-1)/length(clusters),
           x4=nrow(data)*(1:length(clusters))/length(clusters),
           y4=length(genes)+1+labheight*0.1,y3=length(genes)+1+labheight*0.1,y2=length(genes)+1,y1=length(genes)+1)
  poly <- data.frame(id=rep(poly.data$final_id,each=4),color=rep(poly.data$color,each=4))
  poly.x <- numeric()
  poly.y <- numeric()
  for(i in 1:nrow(poly.data)) {
    poly.x <- c(poly.x,poly.data$x1[i],poly.data$x2[i],poly.data$x3[i],poly.data$x4[i])
    poly.y <- c(poly.y,poly.data$y1[i],poly.data$y2[i],poly.data$y3[i],poly.data$y4[i])
  }
  poly <- cbind(poly,poly.x=poly.x,poly.y=poly.y)
  
  # Build the cell type label rectangles
  xlab.rect <- data.frame(xmin=poly.data$x3,
                          xmax=poly.data$x4,
                          ymin=poly.data$y3,
                          ymax=poly.data$y3+labheight*0.9,
                          color=poly.data$color,
                          label=unique(data$final_label))  
  
  # Build the maximum value labels for the right edge
  max.vals <- data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% unlist()
  max.labels <- data.frame(x=nrow(data)+0.01*nrow(data),y=1:length(genes)+0.5,
                           label=sci_label(max.vals))
  max.header <- data.frame(x=nrow(data)+0.075*nrow(data)/2,y=length(genes)+1.5,label="Max data")
  
  # Scale the datas
  data_max <- max(max.vals)
  
  if(logscale) {
    data[genes] <- log10(data[genes] + 1)
    data_max <- log10(data_max)
  }
  
  # Convert to colors
  heat_colors <- colorRampPalette(c("darkblue","dodgerblue","gray80","orangered","red"))(1001)
  
  for(gene in genes) {
    if(normalize_rows == T) {
      data[gene] <- heat_colors[unlist(round(data[gene]/max(data[gene])*1000+1,0))]
    } else {
      data[gene] <- heat_colors[unlist(round(data[gene]/data_max*1000+1,0))]
    }
  }
  
  background_data <- data.frame(xmin=0,xmax=max(data$xpos),ymin=1,ymax=length(genes)+1,fill="#ADCFE0")
  
  # Plot setup
  p <- ggplot(data) +
    scale_fill_identity() +
    theme_classic(base_size=fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    geom_rect(data=background_data,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill))  
  
  # plot the rectangles for each gene
  for(i in 1:length(genes)) {
    
    # plot the rectangles for the heatmap
    p <- p + geom_rect(data=data,aes_string(xmin="xpos-1",xmax="xpos",ymin=i,ymax=i+1,fill=genes[i]))
    
  }
  
  p <- p + 
    # Cluster labels at the top of the plot
    geom_rect(data=xlab.rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=color)) +
    geom_text(data=xlab.rect,aes(x=(xmin+xmax)/2,y=ymin+0.05,label=label),angle=90,vjust=0.35,hjust=0,size=pt2mm(fontsize)) +
    geom_polygon(data=poly,aes(x=poly.x,y=poly.y,fill=color,group=id)) +
    # Maximum value labels at the right edge of the plot
    geom_text(data=max.header,aes(x=x,y=y,label=label),angle=90,hjust=0,vjust=0.5,size=pt2mm(fontsize)) +
    geom_text(data=max.labels,aes(x=x,y=y,label=label),hjust=0,vjust=0.5,size=pt2mm(fontsize))
  
  # Axis scales
  p <- p + scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0),breaks=1:length(genes)+0.5,labels=genes)
  
  return(p)
  
}


#' Violin plots of gene expression for clusters
#' 
#' This function will generate plots similar to Figure 1c of Tasic, et al. (2015).
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param data_source A character object defining where the data is stored. Currently only works with "internal"
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#' 
#' @return a ggplot2 plot object
#' 
#' @examples
#' pottery_plot()
#' 
#' my_genes <- c("Ercc6","Ercc8","Trp53","Pgbd5")
#' my_clusters <- c(1,5,9,10,24,37)
#' pottery_plot(my_genes,my_clusters,logscale=T,fontsize=14)
pottery_plot <- function(genes=c("Hspa8","Snap25","Gad2","Slc17a6"),clusters=1:49,
                         data_source="internal",
                         logscale=F,
                         fontsize=7,labelheight=25) {
  library(dplyr)
  library(ggplot2)
  
  # Fix Ndnf
  genes[genes == "Ndnf"] <- "A930038C07Rik"
  genes[genes == "9630013A20Rik"] <- "X9630013A20Rik"
  
  genes <- rev(genes)

  if(data_source == "internal") {
    data <- scrattch::v1_data
    all.anno <- scrattch::v1_anno
    
    data <- data %>%
      filter(gene %in% genes)
    
    # Reformat the retrieved data
    row.names(data) <- data[,1]
    data <- data %>% 
      select(-1) %>% 
      t() %>% 
      as.data.frame()
    
    data <- data %>%
      mutate(sample_id=row.names(data)) %>%
      select(one_of(c("sample_id",genes)))
    
  }
  
  #Calculate the height of the label:
  labheight <- length(genes)*(labelheight/100)/(1-labelheight/100)
  
  # Build the maximum value labels for the right edge
  max.rect <- data.frame(xmin=length(clusters)+0.5,xmax=length(clusters)+2,
                         ymin=1,ymax=length(genes)+1+labheight)
  max.vals <- data %>% select(-sample_id) %>% summarise_each(funs(max)) %>% unlist()
  max.labels <- data.frame(x=length(clusters)+0.5,y=1:length(genes)+0.5,
                           label=sci_label(max.vals))
  max.header <- data.frame(x=length(clusters)+1.5,y=length(genes)+1,label="Max data")
  
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    filter(final_id %in% clusters) %>%
    arrange(final_id)
  
  # Scale the data
  for(i in 1:length(genes)) {
    gene <- genes[i]
    if(logscale) {
      data[gene] <- log10(data[gene]+1)/log10(max(data[gene])+1)*0.9 + i
    } else {
      data[gene] <- data[gene]/max(data[gene])*0.9 + i
    }  
  }
  
  cluster_positions <- data %>%
    select(final_id) %>%
    unique() %>%
    mutate(xpos=1:n())
  
  data <- data %>% left_join(cluster_positions,by="final_id")
  
  # Variance injection - geom_violin requires some variance, so I add a vanishingly small random number to each data value
  data[genes] <- data[genes] + runif(nrow(data),0,0.00001)
  
  cluster.data <- data %>%
    select(final_label,final_color,final_id,xpos) %>%
    group_by(final_label,final_color,final_id,xpos) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(final_id) %>%
    mutate(labely=length(genes) + 1.1,
           cny=length(genes) + 0.9 + labheight)
  
  hline.frame <- data.frame(y=seq(1,length(genes)+1,1))
  xlab.rect <- data.frame(xmin=seq(0.5,length(clusters)-0.5,1),xmax=seq(1.5,length(clusters)+0.5,1),ymin=length(genes)+1,ymax=length(genes)+1+labheight,color=cluster.data$final_color)
  
  p <- ggplot(data) +
    scale_fill_identity() +
    geom_hline(data=hline.frame,aes(yintercept=y),size=0.2) +
    geom_rect(data=xlab.rect,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=color)) +
    geom_text(data=cluster.data,aes(y=labely,x=xpos,label=final_label),angle=90,hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
    geom_text(data=cluster.data,aes(y=cny,x=xpos,label=cn),size=pt2mm(fontsize)) +
    geom_rect(data=max.rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="white")) +
    geom_text(data=max.header,aes(x=x,y=y,label=label),angle=90,hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
    geom_text(data=max.labels,aes(x=x,y=y,label=label),hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
    scale_y_continuous("",breaks=1:length(genes)+0.5,labels=genes,expand=c(0,0)) +
    scale_x_continuous("",expand=c(0,0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  # plot the violins for each gene
  for(i in 1:length(genes)) {
    p <- p + geom_violin(data=data,aes_string(x="xpos",y=genes[i],fill="final_color"),scale="width",adjust=2)
    p <- p + stat_summary(data=data,aes_string(x="xpos",y=genes[i]),fun.y="median",fun.ymin="median",fun.ymax="median",geom="point",size=1.3)
  }
  
  return(p)
}


#' Boxplots of gene expression for clusters
#' 
#' Generates boxplots for each gene and each cluster. Similar in structure to pottery_plot().
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param data_source A character object defining where the data is stored. Currently only works with "internal"
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#'
#' @return a ggplot2 plot object
#' 
#' @examples
#' boxter_plot()
#' 
#' my_genes <- c("Ercc6","Ercc8","Trp53","Pgbd5")
#' my_clusters <- c(1,5,9,10,24,37)
#' boxter_plot(my_genes,my_clusters,logscale=T,fontsize=14)
boxter_plot <- function(genes=c("Hspa8","Snap25","Gad2","Slc17a6"),clusters=1:49,
                        data_source="internal",
                        logscale=F,
                        fontsize=7,labelheight=25) {
  library(dplyr)
  library(ggplot2)
  
  # Fix Ndnf
  genes[genes == "Ndnf"] <- "A930038C07Rik"  
  genes[genes == "9630013A20Rik"] <- "X9630013A20Rik"
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    data <- scrattch::v1_data
    all.anno <- scrattch::v1_anno
    
    data <- data %>%
      filter(gene %in% genes)
    
    # Reformat the retrieved data
    row.names(data) <- data[,1]
    data <- data %>% 
      select(-1) %>% 
      t() %>% 
      as.data.frame()
    
    data <- data %>%
      mutate(sample_id=row.names(data)) %>%
      select(one_of(c("sample_id",genes)))
    
  }
  
  #Calculate the height of the label:
  labheight <- length(genes)*(labelheight/100)/(1-labelheight/100)
  
  # Build the maximum value labels for the right edge
  max.rect <- data.frame(xmin=length(clusters)+0.5,xmax=length(clusters)+2,
                         ymin=1,ymax=length(genes)+ 1 + labheight)
  max.vals <- data %>% select(-sample_id) %>% summarise_each(funs(max)) %>% unlist()
  max.labels <- data.frame(x=length(clusters)+0.5,y=1:length(genes)+0.5,
                           label=sci_label(max.vals))
  max.header <- data.frame(x=length(clusters)+1.5,y=length(genes)+1,label="Max data")
  
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    filter(final_id %in% clusters) %>%
    arrange(final_id)
  
  # Scale the datas
  for(i in 1:length(genes)) {
    gene <- genes[i]
    if(logscale) {
      data[gene] <- log10(data[gene]+1)/log10(max(data[gene])+1)*0.9 + i
    } else {
      data[gene] <- data[gene]/max(data[gene])*0.9 + i
    }  
  }
  
  cluster_positions <- data %>%
    select(final_id) %>%
    unique() %>%
    mutate(xpos=1:n())
  
  data <- data %>% left_join(cluster_positions,by="final_id")
  
  cluster.data <- data %>%
    select(final_label,final_color,final_id,xpos) %>%
    group_by(final_label,final_color,final_id,xpos) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(final_id) %>%
    mutate(labely=length(genes) + 1.1,
           cny=length(genes) + 0.9 + labheight)
  
  hline.frame <- data.frame(y=seq(1,length(genes)+1,1))
  xlab.rect <- data.frame(xmin=seq(0.5,length(clusters)-0.5,1),xmax=seq(1.5,length(clusters)+0.5,1),ymin=length(genes)+1,ymax=length(genes)+1+labheight,color=cluster.data$final_color)
  
  p <- ggplot(data) +
    scale_fill_identity() +
    geom_hline(data=hline.frame,aes(yintercept=y),size=0.2) +
    geom_rect(data=xlab.rect,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=color)) +
    geom_text(data=cluster.data,aes(y=labely,x=xpos,label=final_label),angle=90,hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
    geom_text(data=cluster.data,aes(y=cny,x=xpos,label=cn,size=40)) +
    geom_rect(data=max.rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="white")) +
    geom_text(data=max.header,aes(x=x,y=y,label=label),angle=90,hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
    geom_text(data=max.labels,aes(x=x,y=y,label=label),hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
    scale_y_continuous("",breaks=1:length(genes)+0.5,labels=genes,expand=c(0,0)) +
    scale_x_continuous("",expand=c(0,0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  # plot the boxplots for each gene
  for(i in 1:length(genes)) {
    p <- p + geom_boxplot(data=data,aes_string(x="xpos",y=genes[i],fill="final_color"))
  }
  
  return(p)
}


#' Heatmaps of gene expression for clusters
#' 
#' Generates heatmaps for each gene and each cluster based on mean gene expression (more functions will be added later).
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param calculation A character object stating which calculation to perform for each heatmap cell. Options: "mean" and "trimmed_mean"
#' @param data_source A character object defining where the data is stored. Currently only works with "internal"
#' @param normalize_rows Logical object, determines if data are normalized to the maximum value for each gene. If FALSE, the heatmap is normalized to the maximum value across all genes.
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#' 
#' @return a ggplot2 plot object
#' 
#' @examples
#' heater_plot()
#' 
#' my_genes <- c("Ercc6","Ercc8","Trp53","Pgbd5")
#' my_clusters <- c(1,5,9,10,24,37)
#' heater_plot(my_genes,my_clusters,logscale=T,fontsize=14)
heater_plot <- function(genes=c("Hspa8","Snap25","Gad2","Slc17a6"),clusters=1:49,
                        calculation="mean",
                        data_source="internal",normalize_rows=FALSE,
                        logscale=T,fontsize=7,toplabel=T,toptext=T,labelheight=25,maxlabel=T,
                        maxval="auto",colorset=c("darkblue","dodgerblue","gray80","orange","orangered")) {
  library(dplyr)
  library(ggplot2)
  
  # Fix Ndnf
  genes[genes == "Ndnf"] <- "A930038C07Rik"  
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    data <- scrattch::v1_data
    all.anno <- scrattch::v1_anno
    
    data <- data %>%
      filter(gene %in% genes)
    
    # Reformat the retrieved data
    row.names(data) <- data[,1]
    data <- data %>% 
      select(-1) %>% 
      t() %>% 
      as.data.frame()
    
    data <- data %>%
      mutate(sample_id=row.names(data)) %>%
      select(one_of(c("sample_id",genes)))
    
    genes[grepl("^[0-9]",genes)] <- paste0("X",genes[grepl("^[0-9]",genes)])
    
    names(data)[2:length(data)] <- genes
  }
  
  #Calculate the height of the label:
  labheight <- length(genes)*(labelheight/100)/(1-labelheight/100)
  
  # Build the maximum value labels for the right edge
  max.rect <- data.frame(xmin=length(clusters)+0.5,xmax=length(clusters)+2,
                         ymin=1,ymax=length(genes) + 1 + labheight)
  max.vals <- data %>% select(-sample_id) %>% summarise_each(funs(max)) %>% unlist()
  max.labels <- data.frame(x=length(clusters)+0.5,y=1:length(genes)+0.5,
                           label=sci_label(max.vals))
  max.header <- data.frame(x=length(clusters)+1.5,y=length(genes)+1.1,label="Max data")
  
  # Filter and order the rows
  data <- left_join(data,all.anno,by="sample_id") %>%
    filter(final_id %in% clusters) %>%
    arrange(final_id)
  
  cluster_positions <- data %>%
    select(final_id) %>%
    unique() %>%
    mutate(xpos=1:n())
  
  data <- data %>% left_join(cluster_positions,by="final_id")
  
  heat_data <- cluster_positions
  
  # Scale the datas
  for(i in 1:length(genes)) {
    gene <- genes[i]
    if(logscale) {
      data[gene] <- log10(data[gene] + 1)
    }
    
    if(calculation == "mean") {
      gene_func <- paste0("mean(",gene,")")
    } else if(calculation == "trimmed_mean") {
      gene_func <- paste0("mean(",gene,",trim = 0.25)")
    } else if(calculation == "percent") {
      gene_func <- paste0("sum(",gene," > 0)/length(",gene,")")
    }
    
    gene_data <- data %>%
      select(one_of("final_id",gene)) %>%
      group_by(final_id) %>%
      summarize_(result = gene_func)
    names(gene_data)[2] <- gene
    
    heat_data <- heat_data %>%
      left_join(gene_data,by="final_id")
  }
  
  # Convert to colors
  heat_colors <- colorRampPalette(colorset)(1001)
  
  if(maxval == "auto") {
    data_max <- max(unlist(heat_data[genes]))
  } else {
    data_max <- maxval
  }
  
  for(gene in genes) {
    if(normalize_rows == T) {
      heat_data[gene] <- heat_colors[unlist(round(heat_data[gene]/max(heat_data[gene])*1000+1,0))]
    } else {
      color_pos <- unlist(round(heat_data[gene]/data_max*1000+1,0))
      color_pos[color_pos > 1001] <- 1001
      heat_data[gene] <- heat_colors[color_pos]
    }
  }
  
  cluster.data <- data %>%
    select(final_label,final_color,final_id,xpos) %>%
    group_by(final_label,final_color,final_id,xpos) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(final_id) %>%
    mutate(labely=length(genes) + 1.1,
           cny=length(genes) + 0.9 + labheight)
  
  xlab.rect <- data.frame(xmin=seq(0.5,length(clusters)-0.5,1),
                          xmax=seq(1.5,length(clusters)+0.5,1),
                          ymin=length(genes) + 1,
                          ymax=length(genes) + 1 + labheight,
                          color=cluster.data$final_color)
  
  p <- ggplot(data) +
    scale_fill_identity() +
    scale_y_continuous("",breaks=1:length(genes)+0.5,labels=genes,expand=c(0,0)) +
    scale_x_continuous("",expand=c(0,0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
    if(toplabel) {
      p <- p +
        geom_rect(data=xlab.rect,aes(xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax,fill=color))
      if(toptext) {
        p <- p +
          geom_text(data=cluster.data,aes(y=labely,x=xpos,label=final_label),angle=90,hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
          geom_text(data=cluster.data,aes(y=cny,x=xpos,label=cn,size=40))
      }
    }
    
    if(maxlabel) {
      p <- p +
        geom_rect(data=max.rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill="white")) +
        geom_text(data=max.header,aes(x=x,y=y,label=label),angle=90,hjust=0,vjust=0.35,size=pt2mm(fontsize)) +
        geom_text(data=max.labels,aes(x=x,y=y,label=label),hjust=0,vjust=0.35,size=pt2mm(fontsize))
    }
    
  
  # plot the heatmap for each gene
  for(i in 1:length(genes)) {
    p <- p + geom_rect(data=heat_data,aes_string(xmin="xpos-0.5",xmax="xpos + 0.5",ymin=i,ymax=i+1,fill=genes[i]))
  }
  
  return(p)
}

testset <- function() eval.parent(substitute( {
  genes <- c("Hspa8","Snap25","Gad2","Slc17a6")
  clusters <- 1:49
  data_source <- "internal"
  sort <- F
  logscale <- F
  fontsize <- 7
} ))