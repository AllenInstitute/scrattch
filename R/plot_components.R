#' Convert font sizes in pt to mm
#' @param pt A numeric font size in pt.
#' @return A numeric font size in mm.
#' 
#' @examples
#' pt2mm(12)
#' 
#' ggplot(mtcars) +
#'   geom_text(aes(x = mpg, y = wt, label = rownames(mtcars)),
#'             size = pt2mm(7))
pt2mm <- function(pt) {
  mm <- pt / 2.834645669
  return(mm)
}

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

#' Split a character string by commas, spaces, and tabs
#' 
#' @param in_string a chr string containing commas, spaces, and/or tabs
#' @return a character vector with each object separated by any combination of commas, spaces, and/or tabs
#' 
#' @examples
#' test <- "Hspa8, Scnn1a,Rbp4    Ptgs2"
#' split_chr(test)
split_cst <- function(in_string) {
  out_chr <- strsplit(in_string,"[, \t]+")[[1]]
  return(out_chr)
}

#' Convert the case of objects in a character vector to Title Case
#' 
#' @param in_chr a character vector
#' @return a character vector with Each Object In Title Case
#' 
#' @examples
#' test <- c("hspa8","scnn1a","fhqwghads")
#' title_case(test)
title_case <- function(in_chr) {
  lower <- tolower(in_chr)
  s <- strsplit(lower, " ")
  result <- paste(toupper(substring(s, 1,1)), substring(s, 2), sep="")
  return(result)
}

#' Convert the case of Riken genes no matter input case
#' 
#' @param in_chr a character vector of Riken gene symbols
#' @return a character vector with correck Riken capitalization
#' 
#' @examples
#' test <- c("6330527o06RiK","A930038C07RIK","a330070k13rik")
#' riken_case(test)
riken_case <- function(in_chr) {
  upper <- toupper(in_chr)
  result <- sub("RIK","Rik",upper)
  return(result)
}

#' Correct case of mouse gene symbols
#' 
#' @param genes a character vector containing gene symbols
#' @return a character vector containing genes with proper case
#' 
#' @examples
#' test <- c("a330070k13rik","HDAC1","sncg","SeRPinB11")
#' fix_mouse_genes(test)
fix_mouse_genes <- function(genes) {
  genes[grepl("[Rr]ik$",genes)] <- riken_case(genes[grepl("[Rr]ik$",genes)])
  genes[!grepl("[Rr]ik$",genes)] <- title_case(genes[!grepl("[Rr]ik$",genes)])
  
  return(genes)
}

#' Convert a character object to a format compatible with SQL SELECT
#' 
#' @param in_chr a character vector
#' @return a character vector consisting of a single string with each object separated by 
#' escaped quotes, bracketted by parens.
#' 
#' @examples
#' test <- c("Hspa8","Cux2","Snap25","Rspo2")
#' chr_to_sql(test)
#' 
#' sql_query <- paste0("SELECT * FROM my_table WHERE gene IN ",chr_to_sql(test))
#' print(sql_query)
chr_to_sql <- function(in_chr) {
  result <- paste0("(",paste0("\"",in_chr,"\"",collapse=","),")")
  return(result)
}

#' Evaluate a character string to a numeric vector
chr_to_num <- function(in_chr) {
  result <- eval(parse(text=paste0("round(c(",in_chr,"),0)")))
  return(result)
}

#' Convert values to colors along a color ramp
#' 
#' @param x a numeric vector to be converted to colors
#' @param minval a number that's used to set the low end of the color scale (default = 0)
#' @param maxval a number that's used to set the high end of the color scale. If NULL (default), 
#' use the highest value in x
#' @param colorset a set of colors to interpolate between using colorRampPalette 
#' (default = c("darkblue","dodgerblue","gray80","orangered","red"))
#' @return a character vector of hex color values generated by colorRampPalette. Color values will
#' remain in the same order as x.
values_to_colors <- function(x, minval = 0, maxval = NULL, colorset = c("darkblue","dodgerblue","gray80","orangered","red")) {

  heat_colors <- colorRampPalette(colorset)(1001)

  if(is.null(maxval)) {
    maxval <- max(x)
  }
  
  heat_positions <- unlist(round((x - minval) / (maxval - minval) * 1000 + 1, 0))
  
  colors <- heat_colors[heat_positions]
  
  colors
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
  
  if(labeltype == "simple") {
    xlab.rect <- data.frame(xmin = 1:nclust - 0.5,
                            xmax = 1:nclust + 0.5,
                            ymin = ngenes + 1,
                            ymax = ngenes + 1 + labheight,
                            color = unique(data$plot_color),
                            label = unique(data$plot_label) )
  }
  
  if(labeltype == "angle") {
    xlab.rect <- data.frame(xmin = (nsamples) * (1:nclust - 1) / nclust,
                            xmax = (nsamples) * (1:nclust) / nclust,
                            # 10% of the label height is reserved for angled polygons
                            ymin = ngenes + 1 + labheight*0.1,
                            ymax = ngenes + 1 + labheight,
                            color = unique(data$plot_color),
                            label = unique(data$plot_label) )
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
