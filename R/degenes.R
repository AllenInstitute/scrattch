#' Find sparse, differentially-expressed genes in pairwise comparisons
#' 
#' This function compares two groups or sets of groups to find genes that are differentially expressed based on
#' the fraction of cells with expression above a given thresshold. Because multiple groups can be combined, these
#' are termed "class 1" and "class 2" for this comparison. The comparison is performed in both directions, so will
#' identify genes with preferrential expression for both class 1 and class 2.
#' 
#' @param data - a data.frame containing the data to be filtered. First column must be "gene", followed by one column per sample.
#' @param anno - a data.frame with annotations for each sample.
#' @param group_by - a character object with the base of the annotation to be used for selecting groups.
#' @param class1_groups - a numeric object with the ids of the groups to be used as class 1
#' @param class2_groups - a numeric object with the ids of the groups to be used as class 2
#' @param data_cut_up - the minimum value cutoff for the "high" expression group. default: 0
#' @param data_cut_dn - the maximum value cutoff for the "low" expression group. default: 0
#' @param frac_cut_up - the minimum fraction of samples in the "high" expression group above data_cut_up. default: 0.95
#' @param frac_cut_dn - the maximum fraction of samples in the "low" expression group above data_cut_dn. default: 0.05
#' @param top - the number of results to keep from each class. default: NULL (returns all results)
#' @param give - either "table" or "genes".
#' 
#' @return Depends on the give parameter: if "table", a data.frame with details for each gene.
#' if "genes", a character vector of the genes that pass the filters.
get_sparse_pairwise_deg <- function(data = v1_data, anno = v1_anno, group_by = "final", 
                                    class1_groups = 1:23, class2_groups = 24:42,
                                    data_cut_up = 0, data_cut_dn = 0,
                                    frac_cut_up = 0.95, frac_cut_dn = 0.05,
                                    top = NULL, give = "genes") {
  library(dplyr)
  
  group_id <- paste0(group_by,"_id")
  
  group1_anno <- anno[anno[,group_id] %in% class1_groups,]
  group1_data <- data %>%
    select(one_of("gene",group1_anno$sample_id))
  
  group2_anno <- anno[anno[,group_id] %in% class2_groups,]
  group2_data <- data %>%
    select(one_of("gene",group2_anno$sample_id))
  
  scores <- data.frame(gene=data$gene)
  
  upcells1 = rowSums(group1_data >  data_cut_up)/ncol(group1_data) #fraction of cells in group 1 "on" for each gene
  upcells2 = rowSums(group2_data >  data_cut_up)/ncol(group2_data) #fraction of cells in group 2 "on" for each gene
  dncells1 = rowSums(group1_data <= data_cut_dn)/ncol(group1_data) #fraction of cells in group 1 "off" for each gene
  dncells2 = rowSums(group2_data <= data_cut_dn)/ncol(group2_data)	#fraction of cells in group 2 "off" for each gene
  
  scores <- cbind(scores,upcells1,upcells2,dncells1,dncells2)
  
  markers1 <- scores %>%
    filter(upcells1 > frac_cut_up & dncells2 >= (1 - frac_cut_dn)) %>%
    arrange(-upcells1)
  
  markers2 <- scores %>%
    filter(upcells2 > frac_cut_up & dncells1 >= (1 - frac_cut_dn)) %>%
    arrange(upcells2)
  
  if(is.null(top)) {
    results <- rbind(markers1,markers2)
  } else {
    results <- rbind(head(markers1,top),tail(markers2,top))
  }
  
  if(give == "genes") {
    results <- as.character(results$gene)
  } else if(give == "table") {
    results <- results
  }
  
  return(results)
  
}