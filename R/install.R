#' Install all packages in the scrattch suite
#'
#' @param packages A character vector of scrattch packages. If NULL (default), installs all packages.
#'
#'
install_scrattch <- function(packages = NULL) {
  if(is.null(packages)) {
    devtools::install_github("AllenInstitute/scrattch.io")
    devtools::install_github("AllenInstitute/iterclust")
    devtools::install_github("AllenInstitute/lowcat")
    devtools::install_github("AllenInstitute/scrattch.vis")
    devtools::install_github("AllenInstitute/tasic2016data")
  } else {
    if("scrattch.io" %in% packages) {
      devtools::install_github("AllenInstitute/scrattch.io")
    }
    if("scrattch.iterclust" %in% packages) {
      devtools::install_github("AllenInstitute/iterclust")
    }
    if("scrattch.lowcat" %in% packages) {
      devtools::install_github("AllenInstitute/lowcat")
    }
    if("scrattch.vis" %in% packages) {
      devtools::install_github("AllenInstitute/scrattch.vis")
    }
    if("tasic2016data" %in% packages) {
      devtools::install_github("AllenInstitute/tasic2016data")
    }
  }
}
