scrattch_dependencies <- list(
  hicat = list(bioconductor = c("impute",
                                "limma",
                                "preprocessCore",
                                "WGCNA"),
               
               github = c("JinmiaoChenLab/Rphenograph")),
  
  io = list(bioconductor = c("rhdf5"),
            
            github = c()),
  
  vis = list(bioconductor = c(),
             
             github = c())
)

#' Install dependencies for scrattch packages from Bioconductor
#'
install_bioc_deps <- function() {
  if(length(scrattch_dependencies$hicat$bioconductor) > 0) {
    BiocManager::install(scrattch_dependencies$hicat$bioconductor)
  }
  
  if(length(scrattch_dependencies$io$bioconductor) > 0) {
    BiocManager::install(scrattch_dependencies$io$bioconductor)
  }
  
  if(length(scrattch_dependencies$vis$bioconductor) > 0) {
    BiocManager::install(scrattch_dependencies$vis$bioconductor)
  }
  
}

#' Install dependencies for scrattch packages from Github
#'
install_github_deps <- function() {
  if(length(scrattch_dependencies$hicat$github) > 0) {
    devtools::install_github(scrattch_dependencies$hicat$github)
  }
  
  if(length(scrattch_dependencies$io$github) > 0) {
    devtools::install_github(scrattch_dependencies$io$github)
  }
  
  if(length(scrattch_dependencies$vis$github) > 0) {
    devtools::install_github(scrattch_dependencies$vis$github)
  }
  
}

#' Install all dependencies from Github and BioConductor
#' 
install_scrattch_deps <- function() {
  install_github_deps()
  install_bioc_deps()
}
