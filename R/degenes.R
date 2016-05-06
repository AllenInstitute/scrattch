#' Find sparse, differentially-expressed genes in pairwise comparisons
#' 
#' @param data - a data.frame containing the data to be filtered. First column must be "gene", followed by one column per sample.
#' @param anno - a data.frame with annotations for each sample.
#' @param group_by - a character object with the base of the annotation to be used for selecting groups.
#' @param class1_group - a numeric object

get_sparse_pairwise_deg <- function(data,anno,group_by,class1_groups,class2_groups,
                                    data_cut_up=0,data_cut_dn=0,
                                    frac_cut_up=0.95,frac_cut_dn=0.05,
                                    top=NULL,give="genes") {
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