#' Get data from GEO
geo_to_list <- function(accession) {
  
  all_packages <- installed.packages()[,"Package"]
  
  if(!sum(c("GEOquery","Biobase") %in% all_packages) == 2) {
    
    run_install <- readline(prompt="Install GEOquery and Biobase from Bioconductor? y/N: \n")
    
    if(run_install == "y") {
      
      install_GEO_deps()
      
    } else {
      stop("Required dependencies not installed. Exiting function.")
    }
  } else {
    
    library(Biobase)
    library(GEOquery)
    library(dplyr)
    
    cat("Retrieving data for",accession,"from GEO\n")
    
    gset <- try(getGEO(accession, GSEMatrix = T), silent = T)
    
    gset <- gset[[1]]
    
    gpl <- annotation(gset)
    
    cat("Retrieving platform annotations for",gpl,"from GEO\n")
    
    platf <- getGEO(gpl, AnnotGPL=TRUE)
    ncbifd <- data.frame(attr(dataTable(platf), "table"))
    
    id_to_names <- ncbifd %>%
      select(ID,Gene.symbol)
    names(id_to_names) <- c("probe_id","gene")
    
    ex <- as.data.frame(exprs(gset))
    
    
    data <- ex %>%
      mutate(probe_id = rownames(ex)) %>%
      left_join(id_to_names) %>%
      select(-probe_id) %>%
      filter(gene != "") %>%
      select(one_of(c("gene",names(ex)))) %>%
      group_by(gene) %>%
      filter(row_number() == 1) %>%
      ungroup() %>%
      mutate(gene = sub("///.+","",gene)) %>%
      mutate(gene = sub("-",".",gene)) %>% 
      arrange(gene) %>%
      as.data.frame()
    
    row.names(data) <- 1:nrow(data)
    
    samples <- names(ex)
    genes <- data$gene
    
    all.anno <- pData(gset)
  
    list(data = data, anno = all.anno)
    
  }
  
}