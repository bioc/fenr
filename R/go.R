#' URL of GO gene association file
#'
#' @return A string with URL.
#' @noRd
get_go_obo_file <- function() {
  getOption("GO_OBO_FILE", "http://purl.obolibrary.org/obo/go.obo")
}

#' URL of GO species webpage
#'
#' @return A string with URL.
#' @noRd
get_go_species_url <- function() {
  getOption("GO_SPECIES_URL", "http://current.geneontology.org/products/pages/downloads.html")
}

#' URL of GO annotation server
#'
#' @return A string with URL.
#' @noRd
get_go_annotation_url <- function() {
  getOption("GO_ANNOTATION_URL", "http://current.geneontology.org/annotations")
}

#' Parse OBO file and return a tibble with key and value
#'
#' @param obo Obo file content as a character vector
#'
#' @return A tibble with term_id, key and value
#' @noRd
parse_obo_file <- function(obo) {
  # Find index of start and end line of each term

  # Start lines
  starts <- stringr::str_which(obo, "\\[Term\\]")

  # Empty lines at end of each term
  blanks <- stringr::str_which(obo, "^$")
  blanks <- blanks[blanks > starts[1]]

  # No space at the end
  if(length(blanks) < length(starts))
    blanks <- c(blanks, length(obo) + 1)

  # End lines: ignore empty lines beyond terms
  ends <- blanks[seq_len(length(starts))]

  # Parse each term
  purrr::map2_dfr(starts, ends, function(i1, i2) {
    obo_term <- obo[(i1 + 1):(i2 - 1)]
    trm <- obo_term |>
      stringr::str_split(":\\s", 2, simplify = TRUE)
    colnames(trm) <- c("key", "value")
    # assuming term_id is in the first line, if not, we are screwed
    tid <- trm[1, 2]
    cbind(trm, term_id = tid) |>
      as.data.frame()
  }) |>
    tibble::as_tibble()
}

#' Extract GO term IDs and names from parsed OBO data
#'
#' @param parsed OBO data parsed by \code{parse_obo_file}
#'
#' @return A tibble with term_id and term_name
#' @noRd
extract_obo_terms <- function(parsed) {
  # Binding variables from non-standard evaluation locally
  key <- term_id <- value <- term_name <- NULL

  terms <- parsed |>
    dplyr::filter(key == "name") |>
    dplyr::select(term_id, term_name = value)

  alt_terms <- parsed |>
    dplyr::filter(key == "alt_id") |>
    dplyr::left_join(terms, by = "term_id") |>
    dplyr::select(term_id = value, term_name)

  dplyr::bind_rows(
    terms,
    alt_terms
  )
}


#' Download GO term descriptions
#'
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `FALSE`.
#'
#' @return A tibble with term_id and term_name.
#' @noRd
fetch_go_terms <- function(use_cache, on_error = "stop") {
  obo_file <- get_go_obo_file()
  if(!assert_url_path(obo_file, on_error))
    return(NULL)

  lpath <- cached_url_path("obo", obo_file, use_cache)
  readr::read_lines(lpath) |>
    parse_obo_file() |>
    extract_obo_terms()
}



#' Find all species available from geneontology.org
#'
#' This function attempts to scrape HTML web page containing a table of
#' available species and corresponding file names. If the structure of the page
#' changes one day and the function stops working, go to
#' \url{http://current.geneontology.org/products/pages/downloads.html} and check
#' file names. The species designation used in this package is the GAF file name
#' without extension (e.g. for a file \file{goa_chicken.gaf} the designation is
#' \file{goa_chicken}).
#'
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns \code{species} and \code{designation}.
#' @export
#' @examples
#' go_species <- fetch_go_species(on_error = "warn")
fetch_go_species <- function(on_error = c("stop", "warn")) {
  on_error <- match.arg(on_error)
  # Binding variables from non-standard evaluation locally
  species <- designation <- `Species/Database` <- File <- NULL

  url <- get_go_species_url()
  qry <- api_query(url, "")
  if(qry$is_error)
    return(catch_error("GO species website", qry$response, on_error))

  u <- qry$response |>
    httr2::resp_body_html() |>
    rvest::html_table()

  u[[1]] |>
    dplyr::mutate(
      species = `Species/Database` |>
        stringr::str_replace_all("\\n", "-") |>
        stringr::str_replace_all("\\s\\s+", " ") |>
        stringr::str_replace_all("(\\S)-", "\\1"),
      designation = File |>
        stringr::str_remove("\\..*$")
    ) |>
    dplyr::select(species, designation)
}


#' Download GO term gene mapping from geneontology.org
#'
#' @param species Species designation. Base file name for species file under
#'   \url{http://current.geneontology.org/annotations}. Examples are
#'   \file{goa_human} for human, \file{mgi} for mouse or \file{sgd} for yeast.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns \code{gene_symbol}, \code{uniprot_id} and \code{term_id}.
#' @noRd
fetch_go_genes_go <- function(species, use_cache, on_error = "stop") {
  # Binding variables from non-standard evaluation locally
  gene_synonym <- db_object_synonym <- symbol <- NULL
  db_id <- go_term <- NULL

  url <- get_go_annotation_url()
  gaf_file <- stringr::str_glue("{url}/{species}.gaf.gz")
  if(!assert_url_path(gaf_file, on_error))
    return(NULL)

  lpath <- cached_url_path(stringr::str_glue("go_gaf_{species}"), gaf_file, use_cache)
  readr::read_tsv(lpath, comment = "!", quote = "", col_names = GAF_COLUMNS,
                  col_types = GAF_TYPES) |>
    dplyr::mutate(gene_synonym = stringr::str_remove(db_object_synonym, "\\|.*$")) |>
    dplyr::select(gene_symbol = symbol, gene_synonym, db_id, term_id = go_term) |>
    dplyr::distinct()
}


#' Get functional term data from gene ontology
#'
#' Download term information (GO term ID and name) and gene-term mapping (gene
#' symbol and GO term ID) from gene ontology.
#'
#' @details This function relies on Gene Ontology's GAF files containing more
#'   generic information than gene symbols. Here, the third column of the GAF
#'   file (DB Object Symbol) is returned as \code{gene_symbol}, but, depending
#'   on the \code{species} argument it can contain other entities, e.g. RNA or
#'   protein complex names. Similarly, the eleventh column of the GAF file (DB
#'   Object Synonym) is returned as \code{gene_synonym}. It is up to the user to
#'   select the appropriate database.
#'
#' @param species Species designation. Examples are \file{goa_human} for human,
#'   \file{mgi} for mouse or \file{sgd} for yeast. Full list of available
#'   species can be obtained using \code{fetch_go_species} - column
#'   \code{designation}.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @importFrom assertthat assert_that
#' @noRd
fetch_go_from_go <- function(species, use_cache, on_error = "stop") {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_species(species, "fetch_go_species", on_error)

  mapping <- fetch_go_genes_go(species = species, use_cache = use_cache, on_error = on_error)
  if(is.null(mapping))
    return(NULL)

  terms <- fetch_go_terms(use_cache = use_cache)

  list(
    terms = terms,
    mapping = mapping
  )
}



#' Download GO term gene mapping from Ensembl
#'
#' @param mart Object class \code{Mart} representing connection to BioMart
#'   database, created with, e.g., \code{useEnsembl}.
#' @param use_cache Logical, if TRUE, the remote data will be cached locally.
#'
#' @return A tibble with columns \code{ensembl_gene_id}, \code{gene_symbol} and
#'   \code{term_id}.
#' @noRd
fetch_go_genes_bm <- function(mart, use_cache = TRUE) {
  # Binding variables from non-standard evaluation locally
  external_gene_name <- term_id <- go_id <- NULL

  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id"),
    mart = mart,
    useCache = use_cache
  ) |>
    dplyr::rename(
      gene_symbol = external_gene_name,
      term_id = go_id
    ) |>
    dplyr::filter(term_id != "") |>
    tibble::as_tibble()
}



#' Get functional term data from Ensembl
#'
#' Download term information (GO term ID and name) and gene-term mapping
#' (gene ID, symbol and GO term ID) from Ensembl.
#'
#' @param mart Object class \code{Mart} representing connection to BioMart
#'   database, created with, e.g., \code{useEnsembl}.
##' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @noRd
fetch_go_from_bm <- function(mart, use_cache = TRUE) {
  assert_that(!missing(mart), msg = "Argument 'mart' is missing.")
  assert_that(is(mart, "Mart"))

  terms <- fetch_go_terms(use_cache = use_cache)
  mapping <- fetch_go_genes_bm(mart, use_cache = use_cache)

  list(
    terms = terms,
    mapping = mapping
  )
}


#' Get Gene Ontology (GO) data
#'
#'
#' This function downloads term information (GO term ID and name) and gene-term
#' mapping (gene ID, symbol, and GO term ID) from either the Ensembl database
#' (using BioMart) or the Gene Ontology database (using GAF files), depending on
#' the provided argument.
#'
#' @details If \code{species} is provided, mapping from a Gene Ontology GAF file
#'   will be downloaded. GAF files contain more generic information than gene
#'   symbols. In this function, the third column of the GAF file (DB Object
#'   Symbol) is returned as \code{gene_symbol}, but, depending on the
#'   \code{species} argument it can contain other entities, e.g. RNA or protein
#'   complex names. Similarly, the eleventh column of the GAF file (DB Object
#'   Synonym) is returned as \code{gene_synonym}. It is up to the user to select
#'   the appropriate database.
#'
#'   Alternatively, if \code{mart} is provided, mapping will be downloaded from
#'   Ensembl database. It will gene symbol and Ensembl gene ID.
#'
#' @param species (Optional) Species designation. Examples are \code{goa_human}
#'   for human, \code{mgi} for mouse, or \code{sgd} for yeast. Full list of
#'   available species can be obtained using \code{fetch_go_species} - column
#'   \code{designation}. This argument is used when fetching data from the Gene
#'   Ontology database.
#' @param mart (Optional) Object class \code{Mart} representing connection to
#'   BioMart database, created with, e.g., \code{useEnsembl}. This argument is
#'   used when fetching data from the Ensembl database.
#' @param use_cache Logical, if TRUE, the remote data will be cached locally.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @importFrom assertthat assert_that
#' @examples
#' # Fetch GO data from Ensembl
#' mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
#' go_data_ensembl <- fetch_go(mart = mart, on_error = "warn")
#' # Fetch GO data from Gene Ontology
#' go_data_go <- fetch_go(species = "sgd", on_error = "warn")
fetch_go <- function(species = NULL, mart = NULL, use_cache = TRUE, on_error = c("stop", "warn")) {
  on_error <- match.arg(on_error)

  assert_that(!(is.null(species) & is.null(mart)),
              msg = "One of the arguments 'species' or 'mart' must be supplied.")
  assert_that(is.null(species) | is.null(mart),
              msg = "Only one of the arguments 'species' or 'mart' must be supplied.")

  if (!is.null(species)) {
    fetch_go_from_go(species, use_cache = use_cache)
  } else  {
    fetch_go_from_bm(mart, use_cache = use_cache)
  }
}
