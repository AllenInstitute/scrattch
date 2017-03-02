annotate_numeric <- function (df, 
                              col = NULL, base = NULL, 
                              scale = "log10", na_val = 0, 
                              colorset = c("darkblue", "white", "red")) 
{
  library(scrattch)
  library(dplyr)
  
  if(is.null(base)) {
    base <- col
  }
  
  if (!is.numeric(df[[col]])) {
    df[[col]] <- as.numeric(df[[col]])
    df[[col]][is.na(df[[col]])] <- na_val
  }
  
  x <- df[[col]]
  
  annotations <- data.frame(label = unique(x)) %>% 
    arrange(label) %>% 
    mutate(id = 1:n())
  
  if (scale == "log10") {
    colors <- values_to_colors(log10(annotations$label + 1), colorset = colorset)
  } else if(scale == "log2") {
    colors <- values_to_colors(log2(annotations$label + 1), colorset = colorset)
  } else if(scale == "zscore") {
    colors <- values_to_colors(scale(annotations$label), colorset = colorset)
  } else {
    colors <- values_to_colors(annotations$label, colorset = colorset)
  }
  annotations <- mutate(annotations, color = colors)
  names(annotations) <- paste0(base, c("_label", "_id", "_color"))
  
  names(df)[names(df) == col] <- paste0(base,"_label")
  df <- left_join(df, annotations, by = paste0(base,"_label"))
  df
}

annotate_categorical <- function(df, 
                                 col = NULL, base = NULL, 
                                 sort_label = T, na_val = "ZZ_Missing", 
                                 colorset = "rainbow", color_order = "sort") {
  
  library(dplyr)
  library(viridis)
  
  if(is.null(base)) {
    base <- col
  }
  
  if(!is.character(df[[col]])) {
    df[[col]] <- as.character(df[[col]])
    df[[col]][is.na(df[[col]])] <- na_val
  }
  
  x <- df[[col]]
  
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
  
  annotations <- mutate(annotations, color = colors)
  
  names(annotations) <- paste0(base,c("_label","_id","_color"))
  names(df)[names(df) == col] <- paste0(base,"_label")
  df <- left_join(df, annotations, by = paste0(base,"_label"))
  df
}

group_annotations <- function(df, keep_order = TRUE) {
  labels <- names(df)[grepl("_label",names(df))]
  if(!keep_order) {
    labels <- labels[order(labels)]
  }
  bases <- sub("_label","",labels)
  
  anno_cols <- c(paste0(rep(bases,each=3),c("_id","_label","_color")))
  extras <- setdiff(names(df),anno_cols)
  
  anno <- select(df,one_of(c("sample_id",anno_cols,extras)))
  
}