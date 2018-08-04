#' Marker genes of various immune cell types.
#'
#' A list of immune cell marker genes compiled from Novershtern et al., Cell (2011).
#'
#' @format A data frame with 15088 rows and 3 variables:
#' \describe{
#'   \item{entrez_id}{Entrez Gene ID}
#'   \item{score}{Cell type score for this gene}
#'   \item{cell_type}{Cell type defined by marker genes}
#' }
#' @source \url{https://www.cell.com/cell/abstract/S0092-8674(11)00005-5}
"novershtern2011"


#' Canonical marker genes of the hematopoietic lineage.
#'
#' A list of canonical immune cell marker genes compiled from Figure 1 in Novershtern et al., Cell (2011).
#'
#' @format A data frame with 135 rows and 5 variables:
#' \describe{
#'   \item{gene_symbol}{Official gene symbol}
#'   \item{entrez_id}{Entrez Gene ID}
#'   \item{score}{Cell type score for this gene}
#'   \item{cell_code}{Coded cell type defined by marker genes}
#'   \item{cell_type}{Named cell type defined by marker genes}
#' }
#' @source \url{https://www.cell.com/cell/abstract/S0092-8674(11)00005-5}
"canonical_immune_markers"

#' Molecular architecture of the mouse nervous system.
#'
#' @format A list with 3 items:
#' \describe{
#'   \item{cell_types}{A data frame describing the 265 nervous cell types.}
#'   \item{cell_markes}{A data frame with marker genes for cell types as symbols and entrez ids.}
#'   \item{cell_metadata}{A data frame with metadata of 160k cells from the mouse nervous tissue.}
#' }
#' @source \url{http://mousebrain.org}
"zeisel2018"
