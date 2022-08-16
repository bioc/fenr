#' Get functional term data from BioPlanet
#'
#' @return A list with terms and mapping tibbles.
#' @export
#'
#' @examples
#' \dontrun{
#' bioplanet_data <- fetch_bp()
#' }
fetch_bp <- function() {
  bp_file <- "https://tripod.nih.gov/bioplanet/download/pathway.csv"
  paths <- readr::read_csv(bp_file, show_col_types = FALSE)

  terms <- paths |>
    dplyr::select(term_id = PATHWAY_ID, term_name = PATHWAY_NAME) |>
    dplyr::distinct()

  mapping <- paths |>
    dplyr::select(term_id = PATHWAY_ID, ncbi_id = GENE_ID, gene_name = GENE_SYMBOL) |>
    dplyr::distinct()

  list(
    terms = terms,
    mapping = mapping
  )
}