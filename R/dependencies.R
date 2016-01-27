#' Retrieve dependencies from Bioconductor
#' 
#' Some packages not available (or not complete) from CRAN can be retrieved from the Bioconductor package repository:
#' WGCNA
#' 
install_bioconductor_dependencies <- function() {
  source("https://bioconductor.org/biocLite.R")
  biocLite("WGCNA")
}