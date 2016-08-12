#' Pull the current version of scrattch from GitHub
update_scrattch <- function() {
  devtools::install_github("hypercompetent/scrattch",auth_token="ed22ec6b1d333fcb9d4d78a1e7ebd29ec72d0048",build_vignettes = T)
}
