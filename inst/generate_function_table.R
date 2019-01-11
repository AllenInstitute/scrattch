library(scrattch)
library(scrattch.io)
library(scrattch.hicat)
library(scrattch.vis)

search_file_for_functions <- function(file,
                                     funs) {
  lines <- readLines(file)
  
  found <- logical(length = length(funs))
  names(found) <- funs
  
  for(i in 1:length(funs)) {
    search_str <- paste0("^",funs[i],"[ <]")
    found[i] <- as.logical(sum(grepl(search_str, lines)))
  }
  
  return(found)
}

make_funs_df <- function(package_name,
                         package_location = NULL) {
  do.call("library", list(package_name))
  
  funs <- lsf.str(paste0("package:",package_name))
  funs_str <- capture.output(print(funs, width = 1e6))
  
  while(sum(grepl(" : ",funs_str)) != length(funs_str)) {
    rm_idx <- integer()
    for(i in 2:length(funs_str)) {
      if(!grepl(" : ", funs_str[i])) {
        funs_str[i - 1] <- paste0(funs_str[i - 1], funs_str[i])
        rm_idx <- c(rm_idx, i)
      }
    }
    funs_str <- funs_str[-rm_idx]
  }
  
  fun_names <- sub(" : .+","",funs_str)
  
  if(!is.null(package_location)) {
    script_dir <- paste0(package_location, "/R")
    if(dir.exists(script_dir)) {
      script_files <- list.files(script_dir, pattern = ".R", full.names = TRUE)
      short_files <- sub(".+/","",script_files)
      
      fun_files <- character(length = length(fun_names))
      
      for(i in 1:length(script_files)) {
        found <- search_file_for_functions(script_files[i], fun_names)
        fun_files[found] <- short_files[i]
      }  
    }
  } else {
    fun_files <- rep("not searched", length(fun_names))
  }
  
  data.frame(package = package_name,
             filename = fun_files,
             name = fun_names,
             parameters = sub(".+function \\(","",sub("\\)[ ]+$","",funs_str)))
}

scrattch_funs <- make_funs_df("scrattch", "/Dropbox/R/scrattch")
scrattch.io_funs <- make_funs_df("scrattch.io", "/Dropbox/R/scrattch.io")
scrattch.hicat_funs <- make_funs_df("scrattch.hicat", "/Dropbox/R/scrattch.hicat")
scrattch.vis_funs <- make_funs_df("scrattch.vis", "/Dropbox/R/scrattch.vis")

all_funs <- rbind(scrattch_funs,
                  scrattch.io_funs,
                  scrattch.hicat_funs,
                  scrattch.vis_funs)
write.csv(all_funs, "inst/scrattch_function_list.csv", row.names = F)
