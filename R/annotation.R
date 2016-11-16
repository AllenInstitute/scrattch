annotate_numeric <- function(x, name = NULL, logscale = F, colorset = c("darkblue","white","red")) {
  
  library(scrattch)
  library(dplyr)
  
  if(!is.numeric(x)) {
    x <- as.numeric(x)
  }
  
  annotations <- data.frame(label = unique(x)) %>%
    arrange(label) %>%
    mutate(id = 1:n())
  
  if(logscale == T) {
    colors <- values_to_colors(log10(annotations$label + 1),colorset = colorset)
  } else {
    colors <- values_to_colors(annotations$label,colorset = colorset)
  }
  
  annotations <- cbind(annotations, color = colors, stringsAsFactors = F)
  
  names(annotations) <- paste0(name,c("_label","_id","_color"))
  
  return(annotations)
  
}

annotate_categorical <- function(x, name = NULL, sort_label = T, colorset = "rainbow", color_order = "sort") {
  
  library(dplyr)
  library(viridis)

  if(!is.character(x)) {
    x <- as.character(x)
  }
  
  annotations <- data.frame(label = unique(x), stringsAsFactors = F)
  
  if(sort_label == T) {
    annotations <- annotations %>% arrange(label)
  }
  
  annotations <- annotations %>%
    mutate(id = 1:n())
  
  if(colorset == "rainbow") {
    colors <- sub("FF$","",rainbow(nrow(annotations)))
  } else if(colorset == "viridis") {
    colors <- sub("FF$","",viridis(nrow(annotations)))
  } else if(colorset == "magma") {
    colors <- sub("FF$","",magma(nrow(annotations)))
  } else if(colorset == "inferno") {
    colors <- sub("FF$","",inferno(nrow(annotations)))
  } else if(colorset == "plasma") {
    colors <- sub("FF$","",plasma(nrow(annotations)))
  } else if(colorset == "terrain") {
    colors <- sub("FF$","",terrain.colors(nrow(annotations)))
  }
  
  if(color_order == "random") {
    
    colors <- sample(colors, length(colors))
    
  }
  
  annotations <- cbind(annotations, color = colors, stringsAsFactors = F)
  
  names(annotations) <- paste0(name,c("_label","_id","_color"))
  
  return(annotations)

}