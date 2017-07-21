#' Read data from a directory of feather files
#' 
get_feather_data <- function(feather_dir, genes, group_by, group_ids) {
  
  library(dplyr)
  library(feather)
  
  data_file <- paste0(feather_dir, "/data.feather")
  anno_file <- paste0(feather_dir, "/anno.feather")
  
  data <- feather(data.file)
  
  # Read annotations and convert factors
  anno <- read_feather(anno_file) %>%
    mutate_if(is.factor, as.character)
  
  # If an _id column was a factor, it's now a character. Convert to numeric for sorting.
  id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != "sample_id"]
  anno[id_cols] <- lapply(anno[id_cols], as.numeric)
  
  # Check the provided genes against the column names in data_file
  data_names <- names(data)
  
  if(sum(genes %in% data_names) != length(genes)) {
    # Report if names don't match after ignorning case
    not_found <- genes[!toupper(genes) %in% toupper(data_names)]
    
    warning(paste(paste0(not_found, collapse = ", "), "not found in feather data!"))
    
    # Update genes to use names as formatted in data
    genes <- data_names[toupper(data_names) %in% toupper(genes)]
  }
  
  # Find column indexes for sample_id and the matched genes
  # This seems to be faster than calling data[,c("sample_id",genes)] directly
  data_cols <- which(data.names %in% c("sample_id", genes))
  
  # Read the data from the data feather file into memory
  gene_data <- data[,feather_cols]
  
  # Change - to . in column names and genes
  colnames(gene_data) <- gsub("-",".",colnames(gene.data))
  genes <- gsub("-",".",genes)
  
  # rename the _id, _label, and _color for the group_by values for use in plotting
  all_anno <- anno %>%
    rename_("plot_id" = paste0(group_by,"_id"),
            "plot_label" = paste0(group_by,"_label"),
            "plot_color" = paste0(group_by,"_color"))
  
  # use the group_ids to retain the order provided by the group_ids argument
  cluster_order <- data.frame(group_ids = group_ids) %>%
    mutate(cluster_x = 1:n())
  
  # Filter and order the rows
  data <- left_join(all_anno, gene.data, by = "sample_id") %>%
    filter(plot_id %in% group_ids) %>%
    left_join(cluster_order, by = c("plot_id" = "group_ids")) %>%
    arrange(cluster_x) %>%
    mutate(xpos = 1:n()) %>%
    select(-plot_id) %>%
    rename_("plot_id" = "cluster_x")
  
  return(data)
}

#' Format data provided in list format for scrattch plots
#' 
#' Currently only compatible with data from feather_to_list()
get_list_data <- function(data_list, genes, group_by, group_ids) {
  
  library(dplyr)
  
  # Read annotations and convert factors
  anno <- data_list$anno %>%
    mutate_if(is.factor, as.character)
  
  # If an _id column was a factor, it's now a character. Convert to numeric for sorting.
  id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != "sample_id"]
  anno[id_cols] <- lapply(anno[id_cols], as.numeric)
  
  # Check the provided genes against the column names in data_file
  data_names <- names(data_list$data)
  
  if(sum(genes %in% data_names) != length(genes)) {
    # Report if names don't match after ignorning case
    not_found <- genes[!toupper(genes) %in% toupper(data_names)]
    
    warning(paste(paste0(not_found, collapse = ", "), "not found in data table!"))
    
    # Update genes to use names as formatted in data
    genes <- data_names[toupper(data_names) %in% toupper(genes)]
  }
  
  gene_data <- data_list$data[,c("sample_id",genes)]
  
  # Change - to . in column names and genes
  colnames(gene_data) <- gsub("-",".",colnames(gene.data))
  genes <- gsub("-",".",genes)
  
  # rename the _id, _label, and _color for the group_by values for use in plotting
  all_anno <- anno %>%
    rename_("plot_id" = paste0(group_by,"_id"),
            "plot_label" = paste0(group_by,"_label"),
            "plot_color" = paste0(group_by,"_color"))
  
  # use the group_ids to retain the order provided by the group_ids argument
  cluster_order <- data.frame(group_ids = group_ids) %>%
    mutate(cluster_x = 1:n())
  
  # Filter and order the rows
  data <- left_join(all_anno, gene.data, by = "sample_id") %>%
    filter(plot_id %in% group_ids) %>%
    left_join(cluster_order, by = c("plot_id" = "group_ids")) %>%
    arrange(cluster_x) %>%
    mutate(xpos = 1:n()) %>%
    select(-plot_id) %>%
    rename_("plot_id" = "cluster_x")
  
  return(data)
}