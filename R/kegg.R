#' KEGG server base URL
#'
#' @return A string with URL. This can be changed by options(KEGG_BASE_URL =
#'   "different/url").
#' @noRd
get_kegg_url <- function() {
  getOption("KEGG_BASE_URL", "https://rest.kegg.jp")
}

#' Find all species available from KEGG
#'
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble, in which column \code{designation} contains species
#'   designations used in function \code{fetch_kegg}.
#' @export
#' @examples
#' spe <- fetch_kegg_species(on_error = "warn")
fetch_kegg_species <- function(on_error = c("stop", "warn")) {
  qry <- api_query(get_kegg_url(), "list/organism")
  if(qry$is_error)
    return(catch_error("KEGG", qry$response, on_error))

  st <- httr2::resp_body_string(qry$response)
  readr::read_tsv(I(st), col_names = c("id", "designation", "species", "phylogeny"), show_col_types = FALSE)
}


#' Download pathway data from KEGG
#'
#' @param species A valid species designation used by KEGG.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A tibble with columns \code{gene_id} and \code{term_id}.
#' @noRd
fetch_kegg_pathways <- function(species, on_error) {
  # Binding variables from non-standard evaluation locally
  term_id <- NULL

  qry <- api_query(get_kegg_url(), stringr::str_glue("list/pathway/{species}"))
  if(qry$is_error)
    return(catch_error("KEGG", qry$response, on_error))

  st <- httr2::resp_body_string(qry$response)
  readr::read_tsv(I(st), col_names = c("term_id", "term_name"), show_col_types = FALSE) |>
    dplyr::mutate(term_id = stringr::str_remove(term_id, "path:"))
}


#' Parse flat file returned from KEGG GET API
#'
#' @param s A string containing flat file result
#'
#' @return A tibble with \code{gene_is}, \code{gene_symbol} and \code{term_id}.
#' @noRd
parse_kegg_genes <- function(s) {
  # Binding variables from non-standard evaluation locally
  data <- gene_symbol <- gene_id <- NULL

  entries <- stringr::str_split(s, "///") |> unlist()
  entries <- entries[-length(entries)] # there is /// at the end
  purrr::map(entries, function(entry) {
    d <- stringr::str_split(entry, "\n") |>
      unlist()
    n <- length(d)
    i <- 1

    # find ENTRY key
    while(i <= n & !(stringr::str_detect(d[i], "^ENTRY")))
      i <- i + 1
    paths <- stringr::str_split(d[i], "\\s+") |> unlist()
    pathway <- paths[2]

    # find GENE key
    i <- 1
    while(i <= n & !(stringr::str_detect(d[i], "^GENE")))
      i <- i + 1
    # no GENE key found
    if(i > n)
      return(tibble::tibble(
        gene_id = character(0),
        gene_symbol = character(0)
      ))

    # extract genes
    genes <- stringr::str_remove(d[i], "GENE\\s+")
    i <- i + 1
    while(i <= n & stringr::str_detect(d[i], "^\\s+")) {
      genes <- c(genes, stringr::str_remove(d[i], "^\\s+"))
      i <- i + 1
    }

    # create final tibble, attempt to extract gene symbols when semicolon is found
    genes |>
      tibble::as_tibble_col(column_name = "data") |>
      tidyr::separate(data, c("gene_id", "gene_symbol"), sep = "\\s+", extra = "merge") |>
      dplyr::mutate(gene_symbol = dplyr::if_else(
        stringr::str_detect(gene_symbol, ";"),
        stringr::str_remove(gene_symbol, ";.+$"),
        gene_id
      )) |>
      tibble::add_column(term_id = pathway)
  }) |>
    purrr::list_rbind()
}


#' Download term - gene mapping from KEGG
#'
#' @param pathways A character vector with KEGG pathways
#' @param batch_size Number of pathways sent to KEGG database in one query. The
#'   maximum allowed is 10.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#' @importFrom assertthat assert_that
#' @return A tibble with columns \code{gene_id} and \code{term_id}
#' @noRd
fetch_kegg_mapping <- function(pathways, batch_size, on_error = "stop") {
  assert_that(is.character(pathways))
  batches <- split(pathways, ceiling(seq_along(pathways) / batch_size))

  raise_error <- FALSE
  pb <- progress::progress_bar$new(total = length(batches))
  tb <- purrr::map(batches, function(batch) {
    dbentries <- paste(batch, collapse = "+")

    # this returns a flat file
    qry <- api_query(get_kegg_url(), stringr::str_glue("get/{dbentries}"))
    if(qry$is_error) {
      catch_error("KEGG", qry$response, on_error)
      raise_error <<- TRUE
      return(NULL)
    }

    pb$tick()
    st <- httr2::resp_body_string(qry$response)
    parse_kegg_genes(st)
  }) |>
    purrr::list_rbind()

  if(raise_error) {
    return(NULL)
  } else {
    return(tb)
  }
}


#' Get functional term data from KEGG
#'
#' Download information (pathway ID and name) and gene-pathway mapping (entrez
#' gene ID, gene symbol and pathway ID) from KEGG. Gene symbols are extracted
#' from gene descriptions. For some species (e.g. yeast), gene symbols are
#' returned instead of entrez IDs and not in gene description.
#'
#' @param species KEGG species code, for example "hsa" for human. The
#'   full list of available KEGG species can be found by using \code{fetch_kegg_species}.
#' @param batch_size Number of pathways sent to KEGG database in one query. The
#'   maximum allowed is 10.
#' @param on_error A character vector specifying the error handling method. It
#'   can take values `"stop"` or `"warn"`. The default is `"stop"`. `"stop"`
#'   will halt the function execution and throw an error, while `"warn"` will
#'   issue a warning and return `NULL`.
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @importFrom assertthat assert_that is.count
#' @export
#' @examples
#' kegg_data <- fetch_kegg("mge", on_error = "warn")
fetch_kegg <- function(species, batch_size = 10, on_error = c("stop", "warn")) {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_that(is.count(batch_size))
  assert_that(batch_size <= 10, msg = "batch_size needs to be between 1 and 10")

  terms <- fetch_kegg_pathways(species, on_error)
  if(is.null(terms))
    return(NULL)
  mapping <- fetch_kegg_mapping(terms$term_id, batch_size, on_error)

  list(
    terms = terms,
    mapping = mapping
  )
}
