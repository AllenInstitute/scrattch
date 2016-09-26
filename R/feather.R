#' Write scrattch data to feather files
build_feather <- function(anno = NULL, data = NULL, desc = NULL, filebase = NULL) {
  
  library(dplyr)
  library(feather)
  
  if(is.null(anno)) {
    stop("anno table required.")
  } else if(is.null(data)) {
    stop("data table required.")
  } else if(is.null(desc)) {
    stop("desc table required.")
  } else {
    
    # Check for old data format
    if(names(data)[1] == "gene") {
    
      # if old format is used, transpose the table
      gn <- data$gene
      data <- t(data[,-1])
      colnames(data) <- gn
      data <- as.data.frame(data)
      data <- cbind(data,sample_id = rownames(data))

    }
    
    all <- left_join(anno, data)
    
    mainfile <- paste0(filebase,"_main.feather")
    cat("Writing main table to",mainfile,"\n")
    write_feather(all, mainfile)
    
    descfile <- paste0(filebase,"_desc.feather")
    cat("Writing desc table to",descfile,"\n")
    write_feather(desc, descfile)
    
  }
  
}

feather_to_list <- function(filebase = NULL, oldformat = F) {
  
  library(feather)
  
  if(is.null(filebase)) {
    stop("filebase required.")
  }
  
  mainfile <- paste0(filebase,"_main.feather")
  descfile <- paste0(filebase,"_desc.feather")
  
  if(!file.exists(mainfile)) {
    stop("main file not found.")
  } else if(!file.exists(descfile)) {
    stop("desc file not found.")
  }
  
  cat("Reading ",descfile,"\n")
  desc <- read_feather(descfile)
  
  cat("Reading ",mainfile,"\n")
  main <- read_feather(mainfile)
  
  cat("Splitting main table to anno and data\n")
  anno_cols <- c("sample_id", paste0(rep(desc$base,each=3),c("_id","_label","_color")))
  data_cols <- names(main)[!names(main) %in% anno_cols]
  
  anno <- main[,anno_cols]
  data <- main[,data_cols]
  
  # Remove any remaining non-numeric columns
  data <- data[,unlist(lapply(data,is.numeric))]
  
  if(oldformat == T) {
    # Transforming for compatibility with get_list_data
    # This step will slow things down substantially.
    # In future, better to change get_list_data.
    cat("Transforming data table for compatibility\n")
    gn <- colnames(data)
    data <- t(data)
    data <- cbind(gene = gn, data)
    colnames(data) <- c("gene",anno$sample_id)
  }
  
  out_list <- list(anno = anno,
                   data = data,
                   desc = desc)
  
  return(out_list)
  
}
