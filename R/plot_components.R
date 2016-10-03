#' Convert integers to scientific notation labels
#' 
#' @param in_num a numeric vector
#' @param sig_figs a number indicating how many significant figures should be displayed.
#' @return a character vector with numeric values reformatted in 1.2E3 format
#' 
#' @examples
#' my_numbers <- c(100,15.359,32687,.000468)
#' 
#' sci_label(my_numbers)
#' 
#' sci_label(my_numbers,sig_figs=3)
sci_label <- function(in_num, sig_figs = 2, type = "plot") {
  labels <- character()
  for(i in 1:length(in_num)) {
    x <- in_num[i]
    if(x == 0) {
      first <- paste0("0", ".", paste0(rep("0", sig_figs - 1), collapse="") )
    } else if(log10(x) %% 1 == 0) {
      first <- substr(x, 1, 1)
      if(sig_figs > 1) {
        first <- paste0(first, ".", paste0(rep("0", sig_figs - 1), collapse=""))
      }
    } else {
      first <- round(x / (10 ^ floor(log10(x))), sig_figs - 1)
    }
    if(x == 0) {
      if(type == "plot") {
        label <- paste0(first, "%*%10^0" )
      } else if(type == "datatable") {
        label <- paste0(first, "\u271510<sup>0</sup>" )
      }
    } else {
      if(type == "plot") {
        label <- paste0(first, "%*%10^", floor(log10(x)))
      } else if(type == "datatable") {
        label <- paste0(first, "\u271510<sup>", floor(log10(x)),"</sup>")
      }
    }
    labels <- c(labels, label)
  }
  return(labels)
}

#' Remove the X-axis (and most other margins)
#' 
#' Makes plots more suitable for use with Illustrator by removing most margins
#' and the X-axis (which is usually replaced by something else in my plots).
#' 
#' To further remove the space below the x-axis, use labs(x = NULL)
#' 
#' Based on theme_classic() from ggplot2.
#' 
#' @examples
#' ggplot(mtcars) +
#'  geom_point(aes(x = mpg, y = wt)) +
#'  theme_no_x_margin() +
#'  labs(x = NULL)
theme_no_x <- function(base_size = 12, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(plot.margin = unit(c(rep(0,4)),"line"),
          axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.margin = unit(0,"cm"),
          axis.ticks.x = element_blank())
}

#' Build polygons from plot data for fancy headers built into the plot area
#' 
build_header_polygons <- function(data, ngenes, nsamples, nclust, labelheight = 25, labeltype = "angle") {
  # Two label types: 
  # "angle" will draw a polygon with the base lined up with samples, and the top
  # divided evenly for each cluster.
  # "square" will draw a rectangle with the top points matching bottom points
  
  library(dplyr)
  
  ## Note on plot dimensions
  # The range of the plot area (not including labels) will be
  # y-axis: 1:(ngenes + 1)
  # x-axis: 0:(nsamples)
  
  # Calculate the height of the label in plot dimensions:
  labheight <- ngenes*(labelheight/100)/(1-labelheight/100)
  
  # Build cell type label polygons
  # polygon points are built in this order: 1 = bottom-right, 2 = bottom-left, 3 = top-left, 4 = top-right
  # For angled labels, the bottom two x positions are calculated based on the number of samples
  # in each cluster. The top positions are evenly spaced based on the number of clusters.
  poly.data <- data %>% 
    group_by(plot_id) %>%
    summarise(color = plot_color[1],
              x1 = max(xpos),
              x2 = min(xpos) - 1) %>%
    mutate(x3 = (nsamples) * (1:nclust - 1) / nclust,
           x4 = (nsamples) * (1:nclust) / nclust,
           # ngenes + 1 is the top of the plot body
           y1 = ngenes + 1,
           y2 = ngenes + 1,
           # The angled portion of the label will be 10% of the total label height 
           y3 = ngenes + 1 + labheight * 0.1,
           y4 = ngenes + 1 + labheight * 0.1)
  
  # For a simpler square label, set the top and bottom x-positions to be the same
  if(labeltype == "square") {
    poly.data <- poly.data %>%
      mutate(x3 = x2,
             x4 = x1)
  }
  
  # Restructure the polygons for ggplot2's geom_poly().
  # The data should have a single x and y column in order, with id and color for each polygon
  poly <- data.frame(id = rep(poly.data$plot_id, each = 4),
                     color = rep(poly.data$color, each = 4))
  poly.x <- numeric()
  poly.y <- numeric()
  for(i in 1:nrow(poly.data)) {
    poly.x <- c(poly.x,poly.data$x1[i],poly.data$x2[i],poly.data$x3[i],poly.data$x4[i])
    poly.y <- c(poly.y,poly.data$y1[i],poly.data$y2[i],poly.data$y3[i],poly.data$y4[i])
  }
  poly <- cbind(poly,poly.x=poly.x,poly.y=poly.y)
  
  poly
}

#' Build colorful, rectangular labels for plot headers in plot space
#' 
build_header_labels <- function(data, ngenes, nsamples, nclust, labelheight = 25, labeltype = "simple") {
  
  # Three label types: 
  # simple, which is for use with cluster-based plots
  # angle, for cell-based plots with "angle"-type polygonal labels
  # square, for cell-based plots with "square"-type labels
  
  ## Note on plot dimensions
  # The range of the plot area (not including labels) will be
  # y-axis: 1:(ngenes + 1)
  # x-axis: 0:(nsamples) (for cell-based plots)
  # x-axis: 1:(nclust + 1) (for cluster-based plots)
  
  labheight <- ngenes*(labelheight/100)/(1-labelheight/100)
  
  data <- data %>%
    select(plot_id,plot_label,plot_color,xpos) %>%
    unique()
  
  if(labeltype == "simple") {
    xlab.rect <- data.frame(xmin = 1:nclust - 0.5,
                            xmax = 1:nclust + 0.5,
                            ymin = ngenes + 1,
                            ymax = ngenes + 1 + labheight,
                            color = data$plot_color,
                            label = data$plot_label )
  }
  
  if(labeltype == "angle") {
    xlab.rect <- data.frame(xmin = (nsamples) * (1:nclust - 1) / nclust,
                            xmax = (nsamples) * (1:nclust) / nclust,
                            # 10% of the label height is reserved for angled polygons
                            ymin = ngenes + 1 + labheight*0.1,
                            ymax = ngenes + 1 + labheight,
                            color = data$plot_color,
                            label = data$plot_label )
  }
  if(labeltype == "square") {
    xlab.rect <- data %>% 
      group_by(plot_id) %>%
      summarise(xmin = min(xpos) - 1,
                xmax = max(xpos),
                ymin = ngenes + 1 + labheight * 0.1,
                ymax = ngenes + 1 + labheight,
                color = plot_color[1],
                label = plot_label[1])
  }

  xlab.rect  
}

#' Covert hclust objects to segments for use in ggplots
#' 
#' @param hc a hclust object
#' @param tree.dir a character object with the direction the tree points to, from root to leaves. options are "down" (default), "up","left", "right".
#' @param dir.lims a 2-member vector with the space in the direction of plotting that the dendrogram will occupy. default = c(0,1)
#' @return a data.frame with segment values for ggplot2's geom_seg. columns: "x","xend","y","yend".
hclust_to_seg <- function(hc, tree.dir = "down", dir.lims = c(0,1)) {
  
  require(ggdendro)
  require(dplyr)
  
  hc.dendro <- as.dendrogram(hc)
  hc.segs <- as.data.frame(segment(dendro_data(hc.dendro)))
  
  ymin = min(dir.lims)
  ymax = max(dir.lims)
  
  yheight = ymax - ymin
  
  norm.segs <- hc.segs %>%
    mutate(y = (y/max(y))*yheight + ymin) %>%
    mutate(yend = (yend/max(yend))*yheight + ymin)
  
  if(tree.dir == "down") {
    
    plot.segs <- norm.segs
    
  } else if(tree.dir == "up") {
    
    ycenter = (ymin + ymax) / 2
    
    plot.segs <- norm.segs %>%
      mutate(y = ycenter + (ycenter - y),
             yend = ycenter + (ycenter - yend))
    
  } else if(tree.dir == "left") {
    
    plot.segs <- norm.segs
    names(plot.segs) <- c("y","x","yend","xend")
    
  } else if(tree.dir == "right") {
    xcenter = (ymin + ymax) / 2
    
    plot.segs <- norm.segs
    names(plot.segs) <- c("y","x","yend","xend")
    
    plot.segs <- plot.segs %>%
      mutate(x = xcenter + (xcenter - x),
             xend = xcenter + (xcenter - xend))
    
  }
  
  plot.segs
  
}

#' Jitter x-y coordinates in a spiral pattern
spiral_jitter <- function(x, y, n = NULL, max_n = NULL, radius = 1, aspect = 1, ratio = "golden") {
  
  if(is.null(n)) {
    n <- length(x)
  }
  
  # pre-calculated golden ratio
  if(ratio == "golden") {
    ratio <- (sqrt(5) + 1) / 2
  }
  
  # calculate the angle based on the ratio
  angle <- 360 / (ratio ^ 2)
  
  # scale the spacing between points
  if(is.null(max_n)) {
    # If no maximum is provided, then c is based
    # on the number of points (n)
    c <- radius / sqrt(n)
  } else {
    # If a maximum is provided, then c is based
    # on the maximum number of points
    c <- radius / sqrt(max_n)
  }
  
  # set up vectors for the jittered results
  x_j <- rep(x,n)
  y_j <- rep(y,n)
  
  # Jitter each point
  for(m in 1:n) {
    # calculate position using polar coordinates
    r <- c * sqrt(m)
    theta <- angle * m
    
    # convert polar coordinates to cartesian coordinates
    if(aspect > 1) {
      # if the aspect is larger than 1 (wider than tall),
      # scale the y positions to compensate
      x_j[m] <- x + r * cos(theta)
      y_j[m] <- y + r * sin(theta) * 1/aspect
    } else {
      # if the aspect is smaller than 1 (taller than wide),
      # scale the x positions to compensate
      x_j[m] <- x + r * cos(theta) * aspect
      y_j[m] <- y + r * sin(theta) 
      
    }
    
  }
  
  results <- data.frame(x = x_j, y = y_j)
  
  return(results)
  
}
