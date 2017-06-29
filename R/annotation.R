#' Generate colors and ids for numeric annotations
#'
#' @param df data frame to annotate
#' @param col name of the numeric column to annotate
#' @param base base name for the annotation, which wil be used in the desc table. If not provided, will use col as base.
#' @param scale The scale to use for assigning colors. Options are "linear","log10","log2, and "zscore"
#' @param na_val The value to use to replace NAs. default = 0.
#' @param colorset A vector of colors to use for the color gradient. default = c("darkblue","white","red")
#' 
#' @return A modified data frame: the annotated column will be renamed base_label, and base_id and base_color columns will be appended
#' 
#' @examples 
#' cars <- mtcars %>%
#'   annotate_num(wt)
#'
#'head(cars)
#'
#'cars2 <- mtcars %>%
#'  annotate_num(wt, weight, "linear", colorset = c("#000000","#808080","#FF0000"))
#'  
#'head(cars2)
annotate_num <- function (df, 
                          col = NULL, base = NULL, 
                          scale = "log10", na_val = 0, 
                          colorset = c("darkblue", "white", "red")) 
{
  library(scrattch)
  library(dplyr)
  
  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }
  
  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }
  
  if (!is.numeric(df[[col]])) {
    df[[col]] <- as.numeric(df[[col]])
  }
  
  df[[col]][is.na(df[[col]])] <- na_val
  
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
  } else if(scale == "linear") {
    colors <- values_to_colors(annotations$label, colorset = colorset)
  }
  annotations <- mutate(annotations, color = colors)
  names(annotations) <- paste0(base, c("_label", "_id", "_color"))
  
  names(df)[names(df) == col] <- paste0(base,"_label")
  df <- left_join(df, annotations, by = paste0(base,"_label"))
  df
}

#' Generate colors and ids for categorical annotations
#'
#' @param df data frame to annotate
#' @param col name of the character column to annotate
#' @param base base name for the annotation, which wil be used in the desc table. If not provided, will use col as base.
#' @param sort_label a logical value to determine if the data in col should be arranged alphanumerically before ids are assigned. default = T.
#' @param na_val The value to use to replace NAs. default = "ZZ_Missing".
#' @param colorset The colorset to use for assigning category colors. Options are "rainbow","viridis","inferno","magma", and "terrain"
#' @param color_order The order in which colors should be assigned. Options are "sort" and "random". "sort" assigns colors in order; "random" will randomly assign colors.
#' 
#' @return A modified data frame: the annotated column will be renamed base_label, and base_id and base_color columns will be appended
#' 
#' @examples 
#' flowers <- iris %>%
#'   annotate_cat(Species)
#'
#'head(flowers)
#'
#'flowers2 <- iris %>%
#'  annotate_num(Species, spp, sort_label = F, colorset = "viridis")
#'  
#'head(flowers2)
annotate_cat <- function(df, 
                         col = NULL, base = NULL, 
                         sort_label = T, na_val = "ZZ_Missing", 
                         colorset = "rainbow", color_order = "sort") {
  
  library(dplyr)
  library(viridis)

  
  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
      stop("Specify a column (col) to annotate.")
  }
  
  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }
  
  if(!is.character(df[[col]])) {
    df[[col]] <- as.character(df[[col]])
  }
  
  df[[col]][is.na(df[[col]])] <- na_val
  
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
  
  names(annotations) <- paste0(base, c("_label","_id","_color"))
  
  names(df)[names(df) == col] <- paste0(base,"_label")
  
  df <- left_join(df, annotations, by = paste0(base, "_label"))
  
  df
}

#'Group annotation columns
#'
#'@param df the annotation dataframe to arrange
#'@param keep_order a logical value. If FALSE, will sort the annotations alphanumerically by base.
#'
#'
#'@return an annotation data frame with reordered columns
#'
#'@examples
#'anno <- mtcars %>%
#'  mutate(sample_id = paste0("car",1:n())) %>%
#'  select(sample_id, wt, mpg) %>%
#'  annotate_num(wt) %>%
#'  annotate_num(mpg)
#'  
#'head(anno)
#'
#'anno2 <- group_annotations(anno)
#'
#'head(anno2)
#'
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