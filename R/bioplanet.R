#' URL of Bioplanet pathway file
#'
#' @return A string with URL.
#' @noRd
get_bioplanet_pathway_file <- function() {
  getOption("BIOPLANET_PATHWAY_FILE", "https://tripod.nih.gov/bioplanet/download/pathway.csv")
}

#' Get functional term data from BioPlanet
#'
#' Download term information (term ID and name) and gene-pathway mapping
#' (NCBI gene ID, gene symbol and pathway ID) from BioPlanet.
#'
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @examples
#' bioplanet_data <- fetch_bp(on_error = "warn")
fetch_bp <- function(use_cache = TRUE, on_error = c("stop", "warn")) {
  on_error <- match.arg(on_error)

  # Binding variables from non-standard evaluation locally
  PATHWAY_ID <- PATHWAY_NAME <- GENE_ID <- GENE_SYMBOL <- NULL

  pathway_file <- get_bioplanet_pathway_file()
  if(!assert_url_path(pathway_file, on_error))
    return(NULL)

  lpath <- cached_url_path("bioplanet_pathway", pathway_file, use_cache)
  paths <- readr::read_csv(lpath, show_col_types = FALSE)

  terms <- paths |>
    dplyr::select(term_id = PATHWAY_ID, term_name = PATHWAY_NAME) |>
    dplyr::distinct()

  mapping <- paths |>
    dplyr::select(term_id = PATHWAY_ID, ncbi_id = GENE_ID, gene_symbol = GENE_SYMBOL) |>
    dplyr::distinct()

  list(
    terms = terms,
    mapping = mapping
  )
}
