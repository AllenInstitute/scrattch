#' Install all packages in the scrattch suite
#'
#' @param packages A character vector of scrattch packages. If NULL (default), installs all packages. Available packages are:
#' \itemize{
#' \item "scrattch.io": File handling for RNA-seq data storage formats.
#' \item "scrattch.hicat": Hierarchical, iterative clustering for analysis of transcriptomics.
#' \item "scrattch.vis": Data visualization for RNA-seq analysis results.
#' \item "tasic2016data": Data and annotations from Tasic, et al., Nature Neuroscience (2016). Used for demonstrations and testing.
#' }
#'
install_scrattch <- function(packages = NULL) {
  install_scrattch_deps()
  
  if(is.null(packages)) {
    devtools::install_github("AllenInstitute/scrattch.io")
    devtools::install_github("AllenInstitute/scrattch.hicat")
    devtools::install_github("AllenInstitute/scrattch.vis")
    devtools::install_github("AllenInstitute/tasic2016data")
  } else {
    if("scrattch.io" %in% packages) {
      devtools::install_github("AllenInstitute/scrattch.io")
    }
    if("scrattch.hicat" %in% packages) {
      devtools::install_github("AllenInstitute/scrattch.hicat")
    }
    if("scrattch.hicat" %in% packages) {
      devtools::install_github("AllenInstitute/scrattch.vis")
    }
    if("tasic2016data" %in% packages) {
      devtools::install_github("AllenInstitute/tasic2016data")
    }
  }
}
