#'Retrieve regions around TSS from the UCSC Genome Browser
#'
#'@param symbols A character vector containing gene symbols
#'@param expand A numeric vector with how many bp to expand upstream and downstream of the TSS. If a single number is supplied, it will be applied both 5' and 3' of the TSS. If 2 numbers are used, the first will be used for 5' and the second for 3'.
#'@param genome A character object with the genome abbreviation e.g. "mm10", "hg19".
#'
#'@return A data frame with BED-format columns chr, start, end, name, score, and strand.
#'
#'@examples
#'
#'my_genes <- c("Hspa8","Sst","Ntsr1","Foxp2")
#'regions_5k <- get_tss_regions(my_genes)
#'
#'regions_1kup_5kdn <- get_tss_regions(my_genes,
#'                                     expand=c(1000,5000),
#'                                     genome="mm9")
get_tss_regions <- function(symbols,expand=5000,genome="mm10") {
  library(dplyr)
  library(rtracklayer)
  
  session <- browserSession("UCSC")
  genome(session) <- genome
  refgene <- getTable(ucscTableQuery(session,table="refGene"))
  
  if(length(expand) == 1) {
    expand5 <- expand
    expand3 <- expand
  } else if(length(expand) == 2) {
    expand5 <- expand[1]
    expand3 <- expand[2]
  }
  
  tss_windows <- refgene %>%
    filter(name2 %in% symbols) %>%
    mutate(tss = ifelse(strand == "+", txStart, txEnd)) %>%
    mutate(tss_up = ifelse(strand == "+", tss - expand5, tss - expand3),
           tss_dn = ifelse(strand == "+", tss + expand3, tss + expand5)) %>%
    select(chrom,tss_up,tss_dn,name2,score,strand)
  
  names(tss_windows) <- c("chr","start","end","name","score","strand")
  
  return(tss_windows)
  
}

#' Convert data.frames in BED-like format to GRanges objects
bed_to_GRanges <- function(bed) {
  library(rtracklayer)
  
  gr <- GRanges(seqnames=bed$chr,
                IRanges(start=bed$start,
                        end=bed$end),
                strand=bed$strand,
                mcols=bed[,c("name","score")])
  
  return(gr)
}