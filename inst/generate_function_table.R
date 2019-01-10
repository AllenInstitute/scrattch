library(scrattch)
library(scrattch.io)
library(scrattch.hicat)
library(scrattch.vis)

make_funs_df <- function(package_name) {
  funs <- lsf.str(paste0("package:",package_name))
  funs_str <- capture.output(print(funs, width = 1e4))
  
  data.frame(package = package_name,
             name = sub(" : .+","",funs_str),
             parameters = sub(".+function \\(","",sub("\\)[ ]+$","",funs_str)))
}

scrattch_funs <- make_funs_df("scrattch")

scrattch.io_funs <- lsf.str("package:scrattch.io")
scrattch.io_funs <- capture.output(print(scrattch.io_funs))

scrattch.hicat_funs <- lsf.str("package:scrattch.hicat")
scrattch.hicat_funs <- capture.output(print(scrattch.hicat_funs))

scrattch.vis_funs <- lsf.str("package:scrattch.vis")
scrattch.vis_funs <- capture.output(print(scrattch.vis_funs))

