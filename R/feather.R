#' Write scrattch data to feather files
build_feather <- function(anno = NULL, data = NULL, desc = NULL, feather_dir = NULL) {
  
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
      data <- as.data.frame(data,stringsAsFactors=F)
      data <- cbind(data,sample_id = rownames(data),stringsAsFactors=F)

    }
    
  if(!dir.exists(feather_dir)) {
    dir.create(feather_dir)
  }
  
  datafile <- paste0(feather_dir,"/data.feather")
  cat("Writing data table to",datafile,"\n")
  write_feather(data, datafile)
  
  annofile <- paste0(feather_dir,"/anno.feather")
  cat("Writing anno table to",annofile,"\n")
  write_feather(anno, annofile)
  
  descfile <- paste0(feather_dir,"/desc.feather")
  cat("Writing desc table to",descfile,"\n")
  write_feather(desc, descfile)
  
  }
}

feather_to_list <- function(feather_dir = NULL, oldformat = F) {
  
  library(feather)
  
  if(is.null(feather_dir)) {
    stop("feather directory required.")
  }
  
  if(!dir.exists(feather_dir)) {
    stop("feather directory doesn't exist")
  }
  
  datafile <- paste0(feather_dir,"/data.feather")
  annofile <- paste0(feather_dir,"/anno.feather")
  descfile <- paste0(feather_dir,"/desc.feather")
  
  if(!file.exists(datafile)) {
    stop("data file not found.")
  } else if(!file.exists(annofile)) {
    stop("anno file not found")
  } else if(!file.exists(descfile)) {
    stop("desc file not found.")
  }
  
  cat("Reading ",descfile,"\n")
  desc <- read_feather(descfile)
  
  cat("Reading ",annofile,"\n")
  desc <- read_feather(annofile)
  
  cat("Reading ",datafile,"\n")
  main <- read_feather(datafile)
  
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
