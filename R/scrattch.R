#' Pull the current version of scrattch from GitHub
update_scrattch <- function() {
  devtools::install_git("http://stash.corp.alleninstitute.org/scm/tran/scrattch.git")
}
