library(scrattch)
library(scrattch.io)
library(scrattch.hicat)
library(scrattch.vis)

make_funs_df <- function(package_name) {
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
  
  data.frame(package = package_name,
             name = sub(" : .+","",funs_str),
             parameters = sub(".+function \\(","",sub("\\)[ ]+$","",funs_str)))
}

scrattch_funs <- make_funs_df("scrattch")
scrattch.io_funs <- make_funs_df("scrattch.io")
scrattch.hicat_funs <- make_funs_df("scrattch.hicat")
scrattch.vis_funs <- make_funs_df("scrattch.vis")

all_funs <- rbind(scrattch_funs,
                  scrattch.io_funs,
                  scrattch.hicat_funs,
                  scrattch.vis_funs)
write.csv(all_funs, "inst/scrattch_function_list.csv", row.names = F)
