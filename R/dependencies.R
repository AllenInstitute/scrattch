#' Install dependencies for scrattch packages from Bioconductor
#'
#' @param packages A character vector of scrattch packages. If NULL (default), installs dependencies for all packages.
#'
#'
install_bioc_deps <- function(packages = NULL) {
  if(is.null(packages)) {
    cat("Installing Bioconductor dependencies for scrattch.io")
    scrattch.io::install_bioc_deps()
    cat("Installing Bioconductor dependencies for scrattch.iterclust")
    scrattch.iterclust::install_bioc_deps()
    cat("Installing Bioconductor dependencies for scrattch.lowcat")
    scrattch.lowcat::install_bioc_deps()
  } else {
    if("scrattch.io" %in% packages) {
      cat("Installing Bioconductor dependencies for scrattch.io")
      scrattch.io::install_bioc_deps()
    }
    if("scrattch.iterclust" %in% packages) {
      cat("Installing Bioconductor dependencies for scrattch.iterclust")
      scrattch.iterclust::install_bioc_deps()
    }
    if("scrattch.lowcat" %in% packages) {
      cat("Installing Bioconductor dependencies for scrattch.lowcat")
      scrattch.lowcat::install_bioc_deps()
    }
  }
}
