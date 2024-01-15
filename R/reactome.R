#' Reactome server base URL
#'
#' @return A string with URL.
#' @noRd
get_reactome_url <- function() {
  getOption("REACTOME_BASE_URL", "https://reactome.org/ContentService")
}

#' URL of Reactome Ensembl file
#'
#' @return A string with URL.
#' @noRd
get_reactome_ensembl_file <- function() {
  getOption("REACTOME_ENSEMBL_FILE", "https://reactome.org/download/current/Ensembl2Reactome.txt")
}

#' URL of Reactome gene association file
#'
#' @return A string with URL.
#' @noRd
get_reactome_gaf_file <- function() {
  getOption("REACTOME_GAF_FILE", "https://reactome.org/download/current/gene_association.reactome.gz")
}

#' List of available Reactome species
#'
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#' @return A tibble with species names used by Reactome.
#' @export
#' @examples
#' re <- fetch_reactome_species(on_error = "warn")
fetch_reactome_species <- function(on_error = c("stop", "warn")) {
  on_error <- match.arg(on_error)

  # Binding variables from non-standard evaluation locally
  dbId <- displayName <- taxId <- NULL

  qry <- api_query(get_reactome_url(), "data/species/main")
  if(qry$is_error)
    return(catch_error("Reactome", qry$response, on_error))

  httr2::resp_body_json(qry$response) |>
    purrr::map(tibble::as_tibble) |>
    purrr::list_rbind() |>
    dplyr::select(db_id = dbId, designation = displayName, tax_id = taxId) |>
    dplyr::distinct()
}


#' Download pathway data from Reactome
#'
#' @param tax_id Taxonomy ID of the species, a string
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns \code{term_id} and \code{term_name}
#' @noRd
fetch_reactome_pathways <- function(tax_id, on_error = "stop") {
  # Binding variables from non-standard evaluation locally
  stId <- displayName <- NULL

  path <- "data/schema/Pathway/min"
  params <- list(
    species = tax_id,
    page = 1,
    offset = 20000
  )

  qry <- api_query(get_reactome_url(), path, params)
  if(qry$is_error)
    return(catch_error("Reactome", qry$response, on_error))

  httr2::resp_body_json(qry$response) |>
    purrr::map(tibble::as_tibble) |>
    purrr::list_rbind() |>
    dplyr::select(term_id = stId, term_name = displayName)
}


#' Download term - Ensembl gene ID mapping from Reactome
#'
#' @details This function downloads one large file containing a mapping between
#'   Enxembl gene IDs and Reactome terms. This is significantly faster than
#'   \code{fetch_reactome_genes}, but if gene symbols are required, needs
#'   additional ID conversion.
#'
#' @param spec Reactome species.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns \code{gene_id} and \code{term_id}
#' @noRd
fetch_reactome_ensembl_genes <- function(spec, use_cache = TRUE, on_error = "stop") {
  # Binding variables from non-standard evaluation locally
  species <- gene_id <- term_id <- NULL

  ensembl_file <- get_reactome_ensembl_file()
  if(!assert_url_path(ensembl_file, on_error))
    return(NULL)

  lpath <- cached_url_path("ensembl2reactome", ensembl_file, use_cache)
  colms <- c("gene_id", "term_id", "url", "event", "evidence", "species")
  readr::read_tsv(lpath, col_names = colms, show_col_types = FALSE) |>
    dplyr::filter(species == spec) |>
    dplyr::select(gene_id, term_id) |>
    dplyr::distinct()
}

#' Download term - gene association from Reactome
#'
#' @param tax_id Taxon ID
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns \code{accession_number}, \code{gene_symbol} and \code{term_id}
#' @noRd
fetch_reactome_gene_association <- function(tax_id, use_cache = TRUE, on_error = "stop") {
  # Binding variables from non-standard evaluation locally
  symbol <- taxon <- db_ref <- db_id <- NULL

  gaf_file <- get_reactome_gaf_file()
  if(!assert_url_path(gaf_file, on_error))
    return(NULL)

  lpath <- cached_url_path("reactome_gaf", gaf_file, use_cache)
  readr::read_tsv(lpath, comment = "!", quote = "", col_names = GAF_COLUMNS, col_types = GAF_TYPES, skip = 4) |>
    dplyr::mutate(
      symbol = stringr::str_remove(symbol, "_.+$"),
      taxon = stringr::str_remove(taxon, "taxon:"),
      db_ref = stringr::str_remove(db_ref, "REACTOME:")
    ) |>
    dplyr::filter(stringr::str_detect(db_ref, "^R-") & taxon == tax_id) |>
    dplyr::select(accession_number = db_id, gene_symbol = symbol, term_id = db_ref) |>
    dplyr::distinct()
}


#' Download term - gene symbol mapping from Reactome
#'
#' @details This function interrogates Reactome API to get term-gene mapping for
#'   all pathways. This is considerable slower than
#'   \code{fetch_reactome_ensembl_genes}. Warning, occasionally, for some
#'   pathways, Reactome does not return gene symbol - only UniProt accession
#'   number is available.
#'
#' @param pathways A character vector with Reactome patway IDs to get
#'   corresponding genes from.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns\code{term_id}, \code{accession_number} and
#'   \code{gene_symbol}.
#' @noRd
fetch_reactome_api_genes <- function(pathways, on_error) {
  # Binding variables from non-standard evaluation locally
  database_name <- gene_symbol <- NULL

  raise_error <- FALSE
  pb <- progress::progress_bar$new(total = length(pathways))
  tb <- purrr::map(pathways, function(pathway) {
    pb$tick()
    path <- stringr::str_glue("data/participants/{pathway}/referenceEntities")
    qry <- api_query(get_reactome_url(), path)
    if(qry$is_error) {
      raise_error <<- TRUE
      return(catch_error("Reactome", qry$response, on_error))
    }

    httr2::resp_body_json(qry$response) |>
      purrr::map(function(js) {
        tibble::tibble(
          database_name = js$databaseName,
          gene_symbol = js$geneName,
          accession_number = js$identifier
        )
      }) |>
      purrr::list_rbind() |>
      dplyr::filter(database_name == "UniProt") |>
      dplyr::select(-database_name) |>
      tidyr::unnest(gene_symbol) |>
      tibble::add_column(term_id = pathway, .before = 1)
  }) |>
    purrr::list_rbind()
  if(raise_error) {
    return(NULL)
  } else {
    return(tb)
  }
}


#' Get functional term data from Reactome
#'
#' Download term information (pathway ID and name) and gene-pathway mapping
#' (Ensembl gene ID or gene symbol and pathway ID) from Reactome.
#'
#' @details Reactome makes mapping between Ensembl ID and pathway ID available
#'   in form of one downloadable file. Also, a gene association file with
#'   mapping between UniProt accession number, gene symbol and Reactome term is
#'   available.  If \code{source = "ensembl"} or \code{source =
#'   "gene_association"} is set, one large file will be downloaded and parsed.
#'   If \code{source = "api"} is set, then Reactome APIs will be interrogated
#'   for each pathway available. This method is considerably slower, especially
#'   for large genomes. However, gene association file contains far fewer
#'   mappings than can be extracted using API.
#'
#' @param species Reactome species designation, for example "Homo sapiens" for
#'   human. Full list of available species can be found using
#'   \code{fetch_reactome_species()}.
#' @param source How to download the mapping. If 'ensembl' or
#'   'gene_association', one mapping file provided by Reactome will be
#'   downloaded, if 'api', then Reactome API will be used. See details.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' reactome_data <- fetch_reactome("Saccharomyces cerevisiae", on_error = "warn")
fetch_reactome <- function(species, source = c("ensembl", "api", "gene_association"),
                           use_cache = TRUE, on_error = c("stop", "warn")) {
  source <- match.arg(source)
  on_error <- match.arg(on_error)
  assert_that(!missing(species), msg = "Argument 'species' is missing.")

  tax_id <- match_species(species, "fetch_reactome_species", "tax_id", on_error)
  if(is.null(tax_id))
    return(NULL)

  terms <- fetch_reactome_pathways(tax_id, on_error)

  if (source == "ensembl") {
    mapping <- fetch_reactome_ensembl_genes(spec = species, use_cache = use_cache)
  } else if (source == "gene_association") {
    mapping <- fetch_reactome_gene_association(tax_id = tax_id, use_cache = use_cache)
  } else {
    mapping <- fetch_reactome_api_genes(terms$term_id, on_error)
  }

  list(
    terms = terms,
    mapping = mapping
  )
}
