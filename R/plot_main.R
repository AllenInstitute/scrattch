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
sample_bar_plot <- function(genes = c("Hspa8","Snap25","Gad2","Vip"),
                         grouping = "final", clusters = 1:10,
                         data_source = "internal",
                         sort = F, logscale = F,
                         fontsize = 7, labelheight = 25, labeltype = "angle",
                         bgcolor = "#ADCFE0") {
  
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    
    # get_internal_data() from data_formatting.R
    data <- get_internal_data(genes,grouping,clusters)
    
  } else if (is.list(data_source)) {
    
    # get_list_data() from data_formatting.R
    data <- get_list_data(data_source,genes,grouping,clusters)
    
  } else {
    
    # get_db_data() from data_formatting.R
    data <- get_db_data(data_source,genes,grouping,clusters)
    
  }
  
  genes <- sub("-",".",genes)
  
  # Calculate the number of genes and clusters for use as plot dimensions
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- data %>% 
    select( one_of(genes) ) %>% 
    summarise_each( funs(max) ) %>% 
    unlist()
  
  # Scale the data based on the logscale option
  # and scale the y-axis values from i to i*0.9 to separate each gene to a row
  for(i in 1:ngenes) {
    gene <- genes[i]
    
    if(logscale) {
      data[gene] <- log10(data[gene] + 1)
    } 
    
    data[gene] <- data[gene] / max(data[gene]) * 0.9 + i
    
  }

  
  # build_header_polygons from plot_components.R
  header_polygons <- build_header_polygons(data = data, ngenes = ngenes, nsamples = nsamples, nclust = nclust, labelheight = labelheight, labeltype = labeltype)
  
  # Build the cell type label rectangles
  header_labels <- build_header_labels(data = data, ngenes = ngenes, nsamples = nsamples, nclust = nclust, labelheight = labelheight, labeltype = labeltype)
  
  # Calculate Plot Scale bars
  scale_bars <- data.frame(gene = genes,
                           ymin = 1:ngenes,
                           ymid = 1:ngenes + 0.45,
                           ymax = 1:ngenes + 0.9,
                           xmin = -nrow(data) * 0.01,
                           xmax = 0)
  
  # Calculate segments for cluster sepration lines
  segment_lines <- data %>%
    group_by(plot_id) %>%
    summarise(x = max(xpos) )

  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = nsamples * 1.01,
                           y = 1:ngenes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = nsamples * 1.01,
                           y = ngenes + 1,
                           label = "Max value")
  

  # The background of the plot is a rectangular object.
  background_data <- data.frame(xmin = 0, 
                                xmax = nsamples, 
                                ymin = 1, 
                                ymax = ngenes + 1, 
                                fill = bgcolor)
  
  # pt2mm function is used for text labels, from plot_components.R
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = 1:ngenes + 0.45, labels = genes) +
    theme_classic(base_size = fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    geom_rect(data = background_data,
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill)) +
    geom_segment(data = segment_lines,
                 aes(x = x, xend = x, y = 1, yend = ngenes + 1),
                 size = 0.2, color = "gray60", linetype = "dashed")
  
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
    geom_rect(data = header_labels,
              aes(xmin = xmin , xmax = xmax, ymin = ymin, ymax = ymax, fill = color) ) +
    geom_text(data = header_labels, 
              aes(x = (xmin + xmax) / 2, y = ymin + 0.05, label = label),
              angle = 90, vjust = 0.35, hjust = 0, size = pt2mm(fontsize)) +
    geom_polygon(data = header_polygons,
                 aes(x = poly.x, y = poly.y, fill = color, group = id) ) +
    # Scale bar elements
    geom_hline(data = scale_bars,
               aes(yintercept = ymin), size = 0.2) +
    geom_segment(data = scale_bars, 
                 aes(x = xmin, xend = xmax,y = ymid, yend = ymid),
                 size = 0.2) +
    geom_segment(data = scale_bars,
                 aes(x = xmin, xend = xmax, y = ymax, yend = ymax),
                 size = 0.2) +
    geom_segment(data = scale_bars,
                 aes(x = xmax, xend = xmax, y = ymin, yend = ymax),
                 size = 0.2) +
    # Maximum value labels at the right edge of the plot
    geom_rect(aes(xmin = nsamples + 1, xmax = (nsamples + 1)*1.15, ymin = 1, ymax = max(header_labels$ymax)), 
              fill = "#FFFFFF") +
    geom_text(data = max_header,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 0, vjust = 0.5, size = pt2mm(fontsize) ) +
    geom_text(data = max_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5, size = pt2mm(fontsize) , parse = TRUE)
  
  p
  
}

#' Heatmaps of gene expression of individual cells
#' 
#' This function will generate plots similar to those in the shiny heatmap geneterator.
#' Warning: this is currently only able to work with internally-supplied datasets (v1_data and v1_anno).
#' Extension to user-supplied datasets will come soon.
#' 
#' @param genes A character vector containing gene symbols to be plotted
#' @param grouping A character object containing the annotation to group data by
#' @param clusters A numeric vector containing clusters to plot (for v1_anno, the range is 1:49)
#' @param data_source A character object defining where the data is stored. Can be a Sqlite3 database file or "internal".
#' @param normalize_rows Logical object, determines if data are normalized to the maximum value for each gene. If FALSE, the heatmap is normalized to the maximum value across all genes.
#' @param logscale Logical object, determines if data is log scaled before plotting.
#' @param fontsize numeric object, the font size (in pts) used to make the plot.
#' @param labelheight numeric object, Percent of the plot height that should be used for the labels (0 to 100).
#' @param labeltype A character object, either "angle" or "square".
#' 
#' @return a ggplot2 plot object
#' 
#' @examples
#' heatcell_plot()
#' 
#' my_genes <- c("Ercc6", "Ercc8", "Trp53", "Pgbd5")
#' my_clusters <- c(1, 5, 9, 10, 24, 37)
#' my_heatcell_plot <- heatcell_plot(my_genes, my_clusters, norm=T, font=12)
#' 
#' ggsave("plot_output.pdf", heatcell_plot, height = 0.2 * length(my_genes) + 2, width = 4)
#' 
#' gene_text <- "Slc17a6 gad2 tac1,RBP4"
#' gene_fix <- fix_mouse_genes(split_cst(gene_text))
#' cluster_text <- "18:22,3,8"
#' cluster_fix <- chr_to_num(cluster_text)
#' 
#' my_heatcell_plot_2 <- heatcell_plot(gene_fix, clust = cluster_fix, font=12)
sample_heatmap_plot <- function(genes = c("Hspa8","Snap25","Gad2","Vip"),
                          clusters = 1:10, grouping = "final",
                          data_source = "internal",
                          logscale = T, normalize_rows = F,
                          fontsize = 7, labelheight = 25,
                          labeltype = "angle") {
  
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    
    data <- get_internal_data(genes,grouping,clusters)
    
  } else {
    
    data <- get_db_data(data_source,genes,grouping,clusters)
    
  }
  
  genes <- sub("-",".",genes)
  
  # Calculate the number of genes and clusters for use as plot dimensions
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals <- data %>% 
    select( one_of(genes) ) %>% 
    summarise_each( funs(max) ) %>% 
    unlist()
  # Calculate the overall maximum value of gene expression
  data_max <- max(max_vals)
  
  # Scale the data and the max based on the logscale option
  if(logscale) {
    data[genes] <- log10(data[genes] + 1)
    data_max <- log10(data_max + 1)
  }
  
  # Convert the data values to heatmap colors
  # values_to_colors() is from plot_components.R
  for(gene in genes) {
    if(normalize_rows == T) {
      data[gene] <- values_to_colors(x = data[gene])
    } else {
      data[gene] <- values_to_colors(x = data[gene], maxval = data_max)
    }
  }
  
  # build_header_polygons from plot_components.R
  header_polygons <- build_header_polygons(data = data, ngenes = ngenes, nsamples = nsamples, nclust = nclust, labelheight = labelheight, labeltype = labeltype)
  
  # Build the cell type label rectangles
  header_labels <- build_header_labels(data = data, ngenes = ngenes, nsamples = nsamples, nclust = nclust, labelheight = labelheight, labeltype = labeltype)

  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = nsamples * 1.01,
                           y = 1:ngenes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = nsamples * 1.01,
                           y = ngenes + 1,
                           label = "Max value")
  
  # Plot setup
  p <- ggplot(data) +
    scale_fill_identity() +
    theme_classic(base_size = fontsize) +
    theme(axis.text = element_text(size=rel(1)),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = 1:ngenes + 0.45, labels = genes)

  # plot the rectangles for each gene
  for(i in 1:length(genes)) {
    
    # plot the rectangles for the heatmap
    p <- p + geom_rect(data = data,
                       aes_string(xmin = "xpos - 1", xmax = "xpos", ymin = i,ymax = i + 1, fill = genes[i]))
    
  }
  
  # Label elements
  # pt2mm() is in plot_components.R
  p <- p + 
    # Cluster labels at the top of the plot
    geom_rect(data = header_labels,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, y = ymin + 0.05, label = label),
              angle = 90, vjust = 0.35, hjust = 0, size = pt2mm(fontsize)) +
    geom_polygon(data = header_polygons,
                 aes(x = poly.x, y = poly.y, fill = color, group = id)) +
    # Maximum value labels at the right edge of the plot
    geom_rect(aes(xmin = nsamples + 1, xmax = (nsamples + 1)*1.15, ymin = 1, ymax = max(header_labels$ymax)), 
              fill = "#FFFFFF") +
    geom_text(data = max_header,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 0, vjust = 0.5, size = pt2mm(fontsize)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.5, size = pt2mm(fontsize), parse = TRUE)
  
  p
  
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
group_violin_plot <- function(genes = c("Hspa8","Snap25","Gad2","Vip"),
                         grouping = "final", clusters = 1:10,
                         data_source = "internal",
                         sort = F, logscale = F,
                         fontsize = 7, labelheight = 25) {
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    
    data <- get_internal_data(genes,grouping,clusters) %>%
      select(-xpos) %>% mutate(xpos = plot_id)
    
  } else {
    
    data <- get_db_data(data_source,genes,grouping,clusters) %>%
      select(-xpos) %>% mutate(xpos = plot_id)
    
  }
  
  genes <- sub("-",".",genes)
  
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  
  # Compute maximum values before scaling to plot space
  max_vals <- data %>% 
    select(one_of(genes)) %>% 
    summarise_each(funs(max)) %>% 
    unlist()
  
  # Variance injection
  # geom_violin() requires some variance, so I add a vanishingly small random number to each data value
  data[genes] <- data[genes] + runif(nrow(data),0,0.00001)
  
  # Scale the data between i and i + 0.9
  for(i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if(logscale) {
      data[gene] <- log10(data[gene] + 1) / log10(gene_max + 1) * 0.9 + i
    } else {
      data[gene] <- data[gene] / gene_max * 0.9 + i
    }
  }
  
  header_labels <- build_header_labels(data = data, ngenes = ngenes, nsamples = 1, nclust = nclust, labelheight = labelheight, labeltype = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = 1:ngenes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = ngenes + 1,
                           label = "Max value")
  
  cluster_data <- data %>%
    group_by(plot_label,plot_color,plot_id) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(plot_id) %>%
    mutate(labely = ngenes + 1.1,
           cny = max(header_labels$ymax) - 0.1,
           xpos = plot_id)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("", breaks = 1:length(genes) + 0.45, labels = genes, expand = c(0, 0)) +
    scale_x_continuous("", expand = c(0, 0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(ngenes)), size = 0.2)
    
  # plot the violins for each gene
  for(i in 1:length(genes)) {
    p <- p + 
      geom_violin(data = data,
                  aes_string(x = "xpos", y = genes[i], fill = "plot_color"),
                  scale = "width", adjust = 2) +
      stat_summary(data = data,
                   aes_string(x = "xpos", y = genes[i]),
                   fun.y = "median", fun.ymin = "median", fun.ymax = "median", geom = "point", size = 0.7)
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, y = ymin + 0.05, label = label),
              angle = 90, vjust = 0.35, hjust = 0, size = pt2mm(fontsize)) +
    # Cluster counts
    geom_text(data = cluster_data,
              aes(y = cny, x = xpos,label = cn),
              size = pt2mm(fontsize)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = nclust + 0.5, xmax = (nclust + 0.5)*1.15, ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 0, vjust = 0.35, size = pt2mm(fontsize)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.35, size = pt2mm(fontsize), parse = TRUE)
    
  p
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
group_box_plot <- function(genes = c("Hspa8","Snap25","Gad2","Vip"),
                         grouping = "final", clusters = 1:10,
                         data_source = "internal",
                         sort = F, logscale = F,
                         fontsize = 7, labelheight = 25) {
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    
    data <- get_internal_data(genes,grouping,clusters) %>%
      select(-xpos) %>% mutate(xpos = plot_id)
    
  } else {
    
    data <- get_db_data(data_source,genes,grouping,clusters) %>%
      select(-xpos) %>% mutate(xpos = plot_id)
    
  }
  
  genes <- sub("-",".",genes)
  
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  
  # Compute maximum values before scaling to plot space
  max_vals <- data %>% 
    select(one_of(genes)) %>% 
    summarise_each(funs(max)) %>% 
    unlist()
  
  # Scale the data between i and i + 0.9
  for(i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if(logscale) {
      data[gene] <- log10(data[gene] + 1) / log10(gene_max + 1) * 0.9 + i
    } else {
      data[gene] <- data[gene] / gene_max * 0.9 + i
    }
  }
  
  header_labels <- build_header_labels(data = data, ngenes = ngenes, nsamples = 1, nclust = nclust, labelheight = labelheight, labeltype = "simple")
  
  # Build the maximum value labels for the right edge
  max_labels <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = 1:ngenes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = ngenes + 1,
                           label = "Max value")
  
  cluster_data <- data %>%
    group_by(plot_label,plot_color,plot_id) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(plot_id) %>%
    mutate(labely = ngenes + 1.1,
           cny = max(header_labels$ymax) - 0.1,
           xpos = plot_id)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("", breaks = 1:length(genes) + 0.45, labels = genes, expand = c(0, 0)) +
    scale_x_continuous("", expand = c(0, 0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    geom_hline(aes(yintercept = 1:(ngenes)), size = 0.2)
  
  # plot the violins for each gene
  for(i in 1:length(genes)) {
    p <- p + 
      geom_boxplot(data = data,
                  aes_string(x = "xpos", y = genes[i], fill = "plot_color"))

  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, y = ymin + 0.05, label = label),
              angle = 90, vjust = 0.35, hjust = 0, size = pt2mm(fontsize)) +
    # Cluster counts
    geom_text(data = cluster_data,
              aes(y = cny, x = xpos,label = cn),
              size = pt2mm(fontsize)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = nclust + 0.5, xmax = (nclust + 0.5)*1.15, ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 0, vjust = 0.35, size = pt2mm(fontsize)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.35, size = pt2mm(fontsize), parse = TRUE)
  
  p
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
group_heatmap_plot <- function(genes=c("Hspa8","Snap25","Gad2","Vip"),clusters=1:10,
                        grouping = "final",calculation="mean",
                        data_source="internal",normalize_rows=FALSE,
                        logscale=T,fontsize=7,labelheight=25,
                        maxval="auto",colorset=c("darkblue","dodgerblue","gray80","orange","orangered")) {
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  if(data_source == "internal") {
    
    data <- get_internal_data(genes,grouping,clusters) %>%
      select(-xpos) %>% mutate(xpos = plot_id)
    
  } else {
    
    data <- get_db_data(data_source,genes,grouping,clusters) %>%
      select(-xpos) %>% mutate(xpos = plot_id)
    
  }
  
  genes <- sub("-",".",genes)
  
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  
  header_labels <- build_header_labels(data = data, ngenes = ngenes, nsamples = 1, nclust = nclust, labelheight = labelheight, labeltype = "simple")
  
  heat_data <- data %>% select(plot_id,xpos)
  
  # Scale the datas
  for(i in 1:length(genes)) {
    gene <- genes[i]

    
    if(calculation == "mean") {
      gene_func <- paste0("mean(",gene,")")
    } else if(calculation == "trimmed_mean") {
      gene_func <- paste0("mean(",gene,",trim = 0.25)")
    } else if(calculation == "percent") {
      gene_func <- paste0("sum(",gene," > 0)/length(",gene,")")
    } else if(calculation == "median") {
      gene_func <- paste0("stats::median(",gene,")")
    }
    
    gene_data <- data %>%
      select(one_of(c("plot_id",gene))) %>%
      group_by(plot_id) %>%
      summarize_(result = gene_func)
    names(gene_data)[2] <- gene
    
    heat_data <- heat_data %>%
      left_join(gene_data,by="plot_id")
    
    
  }
  
  heat_data <- unique(heat_data)
  
  # Build the maximum value labels for the right edge
  
  max_vals <- heat_data %>% 
    select(one_of(genes)) %>% 
    summarise_each(funs(max)) %>% 
    unlist()
  
  max_labels <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = 1:ngenes + 0.5,
                           label = sci_label(max_vals) )
  max_header <- data.frame(x = (nclust + 0.5) * 1.01,
                           y = ngenes + 1,
                           label = "Max value")
  
  if(logscale) {
    heat_data[genes] <- log10(heat_data[genes] + 1)
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
  
  cluster_data <- data %>%
    group_by(plot_label,plot_color,plot_id) %>%
    summarise(cn=n()) %>%
    as.data.frame(stringsAsFactors=F) %>%
    arrange(plot_id) %>%
    mutate(labely = ngenes + 1.1,
           cny = max(header_labels$ymax) - 0.1,
           xpos = plot_id)

  # Plot Setup  
  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("", breaks = 1:length(genes) + 0.5, labels = genes, expand = c(0, 0)) +
    scale_x_continuous("", expand = c(0, 0)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  # plot the heatmap for each gene
  for(i in 1:length(genes)) {
    p <- p + 
      geom_rect(data = heat_data,
                aes_string(xmin = "xpos - 0.5", xmax = "xpos + 0.5", ymin = i, ymax = i + 1, fill = genes[i]))
  }
  
  # Cluster labels
  p <- p +
    geom_rect(data = header_labels, 
              aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = color)) +
    geom_text(data = header_labels,
              aes(x = (xmin + xmax) / 2, y = ymin + 0.05, label = label),
              angle = 90, vjust = 0.35, hjust = 0, size = pt2mm(fontsize)) +
    # Cluster counts
    geom_text(data = cluster_data,
              aes(y = cny, x = xpos,label = cn),
              size = pt2mm(fontsize)) +
    # Maximum value labels on right side of plot
    geom_rect(aes(xmin = nclust + 0.5, xmax = (nclust + 0.5)*1.15, ymin = 1, ymax = max(header_labels$ymax)),
              fill = "white") +
    geom_text(data = max_header,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 0, vjust = 0.35, size = pt2mm(fontsize)) +
    geom_text(data = max_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, vjust = 0.35, size = pt2mm(fontsize), parse = TRUE)
  
  return(p)
}

#' Build Sankey Plots for up to four annotations
group_river_plot <- function(data_source, 
                             group1 = NULL, 
                             group1_filter = NULL,
                             group2 = NULL, 
                             group2_filter = NULL,
                             group3 = NULL, 
                             group3_filter = NULL,
                             group4 = NULL, 
                             group4_filter = NULL) {
  
  library(networkD3)
  library(dplyr)
  
  nodes <- data.frame(name = character(),
                      color = character())
  
  links <- data.frame(source = numeric(),
                      target = numeric(),
                      value = numeric())
  
  anno <- db_to_list(data_source,get_tables = "anno")$anno
  
  # Filtering annotations
  if(!is.null(group1_filter)) {
    filt <- paste0(group1,"_id %in% ",group1_filter)
    anno <- anno %>%
      filter_(filt)
  }
  if(!is.null(group2_filter)) {
    filt <- paste0(group2,"_id %in% ",group2_filter)
    anno <- anno %>%
      filter_(filt)
  }
  if(!is.null(group3_filter)) {
    filt <- paste0(group3,"_id %in% ",group3_filter)
    anno <- anno %>%
      filter_(filt)
  }
  if(!is.null(group4_filter)) {
    filt <- paste0(group4,"_id %in% ",group4_filter)
    anno <- anno %>%
      filter_(filt)
  }
  
  # Add nodes to nodes table
  if(!is.null(group1)) {
    base <- group1
    anno_id <- paste0(base,"_id")
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    group_nodes <- anno %>%
      select(one_of(anno_id,anno_label,anno_color)) %>%
      unique() %>%
      arrange_(anno_id) %>%
      select(-one_of(anno_id)) %>%
      mutate(group = base)
    
    names(group_nodes) <- c("name","color","group")
    
    nodes <- rbind(nodes,group_nodes)
  }
  
  if(!is.null(group2)) {
    base <- group2
    anno_id <- paste0(base,"_id")
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    group_nodes <- anno %>%
      select(one_of(anno_id,anno_label,anno_color)) %>%
      unique() %>%
      arrange_(anno_id) %>%
      select(-one_of(anno_id)) %>%
      mutate(group = base)
    
    names(group_nodes) <- c("name","color","group")
    
    nodes <- rbind(nodes,group_nodes)
  }
  
  if(!is.null(group3)) {
    base <- group3
    anno_id <- paste0(base,"_id")
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    group_nodes <- anno %>%
      select(one_of(anno_id,anno_label,anno_color)) %>%
      unique() %>%
      arrange_(anno_id) %>%
      select(-one_of(anno_id)) %>%
      mutate(group = base)
    
    names(group_nodes) <- c("name","color","group")
    
    nodes <- rbind(nodes,group_nodes)
  }
  
  if(!is.null(group4)) {
    base <- group4
    anno_id <- paste0(base,"_id")
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    group_nodes <- anno %>%
      select(one_of(anno_id,anno_label,anno_color)) %>%
      unique() %>%
      arrange_(anno_id) %>%
      select(-one_of(anno_id)) %>%
      mutate(group = base)
    
    names(group_nodes) <- c("name","color","group")
    
    nodes <- rbind(nodes,group_nodes)
  }
  
  # Add ID column
  nodes <- nodes %>%
    mutate(id = 1:nrow(nodes) - 1)
  
  # Add Edges to the edges table
  if(!is.null(group1) & !is.null(group2)) {
    source_base <- group1
    target_base <- group2
    
    source_label <- paste0(source_base,"_label")
    target_label <- paste0(target_base,"_label")
    
    group_links <- anno %>%
      group_by_(source_label,target_label) %>%
      summarise(value = n())
    names(group_links) <- c("source_lab","target_lab","value")
    group_links <- group_links %>%
      rowwise() %>%
      mutate(source = nodes$id[nodes$group == source_base & nodes$name == source_lab],
             target = nodes$id[nodes$group == target_base & nodes$name == target_lab]) %>%
      select(source,target,value)
    
    links <- rbind(links, group_links)
  }
  if(!is.null(group2) & !is.null(group3)) {
    source_base <- group2
    target_base <- group3
    
    source_label <- paste0(source_base,"_label")
    target_label <- paste0(target_base,"_label")
    
    group_links <- anno %>%
      group_by_(source_label,target_label) %>%
      summarise(value = n())
    names(group_links) <- c("source_lab","target_lab","value")
    group_links <- group_links %>%
      rowwise() %>%
      mutate(source = nodes$id[nodes$group == source_base & nodes$name == source_lab],
             target = nodes$id[nodes$group == target_base & nodes$name == target_lab]) %>%
      select(source,target,value)
    
    links <- rbind(links, group_links)
  }
  if(!is.null(group3) & !is.null(group4)) {
    source_base <- group3
    target_base <- group4
    
    source_label <- paste0(source_base,"_label")
    target_label <- paste0(target_base,"_label")
    
    group_links <- anno %>%
      group_by_(source_label,target_label) %>%
      summarise(value = n())
    names(group_links) <- c("source_lab","target_lab","value")
    group_links <- group_links %>%
      rowwise() %>%
      mutate(source = nodes$id[nodes$group == source_base & nodes$name == source_lab],
             target = nodes$id[nodes$group == target_base & nodes$name == target_lab]) %>%
      select(source,target,value)
    
    links <- rbind(links, group_links)
  }
  
  # build JavaScript object for colors
  
  d3.colors <- paste0("d3.scale.ordinal().domain([\"",paste(nodes$id,collapse="\",\""),"\"]).range([\"",paste(nodes$color,collapse="\",\""),"\"]);")
  
  # Return the plot
  sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
                Target = 'target', Value = 'value', NodeID = 'name',
                NodeGroup = 'name', colourScale = JS(d3.colors),
                fontSize = 12)
  
}



### Legacy function names
barcell_plot <- sample_bar_plot
heatcell_plot <- sample_heatmap_plot
pottery_plot <- group_violin_plot
heater_plot <- group_heatmap_plot
boxter_plot <- group_box_plot


testset <- function() eval.parent(substitute( {
  genes <- c("Hspa8","Snap25","Gad2","Slc17a6")
  clusters <- 1:49
  data_source <- "internal"
  sort <- F
  logscale <- F
  fontsize <- 7
} ))