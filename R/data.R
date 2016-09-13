#' Annotations for V1 single-cell transcriptomic datasets
#' 
#' A dataset consisting of descriptive annotation data for a subset of the
#' cells characterized in Tasic, et al. (2016).
#' 
#' @format A Data Frame with 11 columns and 282 rows
#' \describe{
#'  \item{sample_id}{A unique cell identifier that uniquely corresponds to a single cell. Matches column names in v1_data}
#'  \item{short_name}{A unique short cell identifier that is used in the Allen Brain Institute Vignette}
#'  \item{cre_id}{Numeric ID corresponding to the Cre line that the cell was harvested from}
#'  \item{cre_color}{Character value for a color corresponding to the cre_id}
#'  \item{cre_label}{Character value for a descriptive label for each cre_id}
#'  \item{final_id}{Numeric ID corresponding to the final cluster assignment for each cell}
#'  \item{final_color}{Character value for a color corresponding to final_id}
#'  \item{final_label}{Character value for a descriptive label for each cluster}
#'  \item{depth_id}{Numeric ID corresponding to the sequencing depth for each cell}
#'  \item{depth_color}{Character value for a color corresponding to sequencing depth}
#'  \item{depth_label}{Character value containing the sequencing depth in millions of reads.}
#' }
"v1_anno"

#' RPKM data for single-cell transcriptomic datasets
#' 
#' A dataset consisting of RPKM values for a subset of the cells
#' characterized in Tasic, et al. (2016).
#' 
#' Each row corresponds to a single gene symbol. Each column corresponds to a single cell.
#' 
#' @format A Data Frame with 283 columns (cells) and 24,057 rows (genes)
#' \describe{
#'  \item{gene}{The first column, a list of gene symbols}
#'  \item{other columns}{Each subsequent column is named using a unique cell identifier that corresponds to the filename column in the v1_anno annotations table.}
#' }
#' 
"v1_data"