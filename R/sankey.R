library(ggplot2)
library(dplyr)
library(feather)

# Error function
# used to generate sigmoidal curve
erf <- function(x) {
  2 * pnorm(x * sqrt(2)) - 1
}

# Generates x-y coordinates for a sigmoidal line from x,y to xend,yend with
# the given number of steps.
# additional arguments define what curve to use
# and how far in the x-direction to use it
sigline <- function(x = 0, xend = 1, 
                    y = 0, yend = 1, 
                    steps = 50,
                    sigfun = "erf",
                    sigx = 1.5) {
  
  xsteps <- seq(-sigx, sigx, length.out = steps)
  if(sigfun == "erf") {
    ysteps <- erf(xsteps)
  }
  
  xsteps <- (xsteps + sigx)/(2*sigx)
  ysteps <- (ysteps + max(ysteps))/(2*max(ysteps))
  
  xscaled <- xsteps*(xend - x) + x
  yscaled <- ysteps*(yend - y) + y
  
  data.frame(x = c(x, xscaled, xend),
             y = c(y, yscaled, yend))
  
}

# expands a sigline into a ribbon by adding a height
# can expand using original line as the top, bottom, or mid-point 
# of the ribbon
sigribbon <- function(sigline, height, from = "top") {
  library(dplyr)
  
  if(from == "top") {
    ribbon <- sigline %>%
      mutate(ymin = y - height)
  } else if(from == "bot") {
    ribbon <- sigline %>%
      rename(y = ymin) %>%
      mutate(y = ymin + height)
  } else if(from == "mid") {
    ribbon <- sigline %>%
      mutate(y = y + height/2,
             ymin = y - height)
  }
  
  ribbon
  
}

make_group_nodes <- function(anno,
                        group_by,
                        xpos = NULL) {
  
  nodes <- data.frame(id = numeric(),
                      name = character(),
                      color = character(),
                      n = numeric(),
                      group = character(),
                      xpos = numeric())
  
  for(i in 1:length(group_by)) {
    base <- group_by[i]
    anno_id <- paste0(base,"_id")
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    grouping <- c(anno_id,anno_label,anno_color)
    
    group_nodes <- anno %>%
      select(one_of(grouping)) %>%
      group_by_(.dots = grouping) %>%
      summarise(n = n()) %>%
      arrange_(anno_id) %>%
      mutate(group = base) %>%
      ungroup()
    
    if(is.null(xpos)) {
      group_nodes <- mutate(group_nodes, xpos = i)
    } else {
      group_nodes <- mutate(group_nodes, xpos = xpos[i])
    }
    
    names(group_nodes) <- c("id","name","color","n","group","xpos")
    
    nodes <- rbind(nodes, group_nodes)
  }
  
  nodes
}

make_plot_nodes <- function(group_nodes,
                            # % of height to distribute for padding between nodes
                            pad = 0.1,
                            # plot space for total width
                            width = 0.1) {
  
  total_n <- sum(group_nodes$n)/length(unique(group_nodes$group))
  total_pad <- total_n * pad
  
  group_rects <- group_nodes %>%
    group_by(group) %>%
    mutate(xmin = xpos - width/2,
           xmax = xpos + width/2,
           n_groups = n(),
           group_pad = ifelse(n_groups > 1, 
                              total_pad/(n() - 1), 
                              total_pad),
           n_cum = cumsum(n),
           ymin = ifelse(n_groups > 1,
                         lag(n_cum, default = 0) + (1:n() - 1)*group_pad,
                         group_pad / 2),
           ymax = ifelse(n_groups > 1,
                         n_cum + (0:(n()-1))*group_pad,
                         group_pad /2 + n))
  
  group_rects
  
}

make_group_links <- function(anno,
                        group_by,
                        plot_nodes) {
  
  pairs <- list()
  
  for(i in 2:length(group_by)) {
    pair <- group_by[(i-1):i]
    pairs <- c(pairs, list(pair))
  }
  
  for(pair in pairs) {
    base <- pair
    anno_id <- paste0(base,"_id")
    anno_label <- paste0(base,"_label")
    anno_color <- paste0(base,"_color")
    
    group1_nodes <- plot_nodes %>%
      filter(group == base[1]) %>%
      select(group, id, xmax, ymin)
    names(group1_nodes) <- c("group1",anno_id[1],"x","group1_min")
    
    group2_nodes <- plot_nodes %>%
      filter(group == base[2]) %>%
      select(group, id, xmin, ymin)
    names(group2_nodes) <- c("group2",anno_id[2],"xend","group2_min")
    
    grouping <- c(anno_id, anno_label, anno_color)
    
    group_links <- anno %>%
      select(one_of(grouping)) %>%
      group_by_(.dots = grouping) %>%
      summarise(n = n()) %>%
      arrange_(.dots = anno_id) %>%
      mutate(group1 = base[1],
             group2 = base[2]) %>%
      ungroup() %>%
      left_join(group1_nodes) %>%
      left_join(group2_nodes) %>%
      group_by_(.dots = anno_id[1]) %>%
      arrange_(.dots = anno_id[2]) %>%
      mutate(y = group1_min + cumsum(n)) %>%
      ungroup() %>%
      group_by_(.dots = anno_id[2]) %>%
      arrange_(.dots = anno_id[1]) %>%
      mutate(yend = group2_min + cumsum(n)) %>%
      ungroup()
      
    names(group_links) <- c("group1_id","group2_id",
                            "group1_label","group2_label",
                            "group1_color","group2_color",
                            "n","group1","group2",
                            "x","group1_min","xend","group2_min",
                            "y","yend")
    
    group_links <- group_links %>%
      rowwise() %>%
      mutate(link_id = paste0(group1_label,"_",group1_id,"_to_",
                              group2_label,"_",group2_id)) %>%
      ungroup()
    
    if(exists("all_links")) {
      all_links <- rbind(all_links, group_links)
    } else {
      all_links <- group_links
    }
    
  }
  
  all_links
  
}

make_plot_links <- function(group_links,
                            fill = NULL) {
  
  for(i in 1:nrow(group_links)) {
    
    link_line <- sigline(x = group_links$x[i], xend = group_links$xend[i],
                         y = group_links$y[i], yend = group_links$yend[i])
    
    link_ribbon <- sigribbon(link_line, h = group_links$n[i])
    
    if(is.null(fill)) {
      link_ribbon <- mutate(link_ribbon, fill = "#808080")
    } else if(fill == group_links$group1[i]) {
      link_ribbon <- mutate(link_ribbon, fill = group_links$group1_color[i])
    } else if(fill == group_links$group2[i]) {
      link_ribbon <- mutate(link_ribbon, fill = group_links$group2_color[i])
    } else {
      link_ribbon <- mutate(link_ribbon, fill = "#808080")
    }
    
    if(exists("all_ribbons")) {
      all_ribbons <- rbind(all_ribbons, link_ribbon)
    } else {
      all_ribbons <- link_ribbon
    }
  }
  
  all_ribbons
  
}

build_river_plot <- function(anno, group_by, pad = 0.1, fill_group = NULL) {
  
  group_nodes <- make_group_nodes(anno, group_by)
  plot_nodes <- make_plot_nodes(group_nodes, pad = pad)
  
  group_links <- make_group_links(anno, group_by, plot_nodes)
  plot_links <- make_plot_links(group_links, fill = fill_group)
  
  p <- ggplot() +
    geom_rect(data = plot_nodes,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = color),
              color = "#808080") +
    geom_ribbon(data = plot_links,
                aes(x = x, ymax = y,
                    ymin = ymin,
                    group = link_id,
                    fill = fill),
                color = "#808080",
                alpha = 0.4) +
    scale_fill_identity() +
    scale_y_reverse() +
    theme_void()
  
  p
  
}

build_river_plot_bokeh <- function(anno, group_by, pad = 0.1, fill_group = NULL) {
  library(rbokeh)
  
  group_nodes <- make_group_nodes(anno, group_by)
  plot_nodes <- make_plot_nodes(group_nodes, pad = pad)
  
  #must be sorted by link_id for polygon hover to work
  group_links <- make_group_links(anno, group_by, plot_nodes) %>%
    arrange(link_id)
  
  poly_links <- data.frame(x = numeric, y = numeric, link_id = character())
  for(i in 1:nrow(group_links)) {
    plot_links <- make_plot_links(group_links[i,], fill = fill_group)
    
    poly_link <- data.frame(x = c(rev(plot_links$x),plot_links$x),
                            y = c(rev(plot_links$ymin),plot_links$y),
                            link_id = rep(plot_links$link_id, 2))
    
    poly_links <- rbind(poly_links, poly_link)
      
  }
  
  # hover behavior is strange - it doesn't follow grouping,
  # so these have to be added in sequence
  poly_links$group1_label <- group_links$group1_label
  poly_links$group2_label <- group_links$group2_label
  poly_links$n <- group_links$n
  
  b <- figure(height = 750, width = 1000) %>%
    ly_rect(data = plot_nodes,
            xleft = xmin, xright = xmax,
            ybottom = -ymin, ytop = -ymax,
            color = color,
            hover = list("N Cells" = n),
            fill_alpha = 1) %>%
    ly_polygons(data = poly_links,
                xs = x, ys = -y,
                group = link_id,
                hover = list("Group 1" = group1_label,
                             "Group 2" = group2_label,
                             "N Cells" = n),
                color = "#808080")
  
  b
  
}

build_river_plot_bokeh(anno, c("region","cre","layer"), pad = 0.1)

anno <- read_feather("//AIBSData/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170405/anno.feather")

filtered_anno <- anno %>%
  filter(cre_label == "Snap25-IRES2-Cre") %>%
  filter(region_label == "VISp")

test <- build_river_plot(filtered_anno, c("cre","layer","cluster"), pad = 0.1)

ggsave("~/../Desktop/river_test.pdf",test, width = 12, height = 8)
  