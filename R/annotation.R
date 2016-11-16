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
  
  results <- data.frame(label = x) %>%
    left_join(annotations) %>%
    select(id,label,color)
  
  names(results) <- paste0(name,c("_id","_label","_color"))
  
  return(results)
  
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
    
    colors <- rainbow(nrow(annotations))
    
  } else if(colorset == "viridis") {
    
    colors <- sub("FF$","",viridis(nrow(annotations)))
    
  }
  
  if(color_order == "random") {
    
    colors <- sample(colors, length(colors))
    
  }
  
  annotations <- cbind(annotations, color = colors, stringsAsFactors = F)
  
  results <- data.frame(label = x, stringsAsFactors = F) %>%
    left_join(annotations) %>% 
    select(id,label,color)
  
  names(results) <- paste0(name,c("_id","_label","_color"))
  
  return(results)
}