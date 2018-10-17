#' Install dependencies for scrattch packages from Bioconductor
#'
install_bioc_deps <- function() {
  source("https://bioconductor.org/biocLite.R")
  
  # Required by scrattch.io: rhdf5
  biocLite("rhdf5")
  
  # Required by scrattch.hicat: impute, limma, WGCNA
  # impute and limma are WGCNA dependencies
  biocLite("impute")
  biocLite("limma")
  biocLite("preprocessCore")
  biocLite("WGCNA")
  
}

#' Install dependencies for scrattch packages from Github
#'
install_github_deps <- function() {
  # Required by scrattch.hicat: JinmiaoChenLab/Rphenograph
  devtools::install_github("JinmiaoChenLab/Rphenograph")
  
}

#' Install all dependencies from Github and BioConductor
#' 
install_scrattch_deps <- function() {
  install_github_deps()
  install_bioc_deps()
}