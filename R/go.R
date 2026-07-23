#' URL of GO obo top dir
#'
#' @return A string with URL.
#' @noRd
get_go_obo_dir <- function() {
  getOption("GO_OBO_DIR", "http://purl.obolibrary.org/obo")
}

#' URL of GO gene association file
#'
#' @return A string with URL.
#' @noRd
get_go_obo_file <- function() {
  getOption("GO_OBO_FILE", "http://purl.obolibrary.org/obo/go.obo")
}

#' URL of GO species metadata file
#'
#' @return A string with URL.
#' @noRd
get_go_species_url <- function() {
  getOption("GO_SPECIES_URL", "https://current.geneontology.org/metadata/goex.yaml")
}

#' URL of GO annotation server
#'
#' @return A string with URL.
#' @noRd
get_go_annotation_url <- function() {
  getOption("GO_ANNOTATION_URL", "http://current.geneontology.org/annotations")
}

#' URL of Ensembl's biomart
#'
#' @return A string with URL.
#' @noRd
get_biomart_url <- function() {
  getOption("ENSEMBL_BIOMART", "http://www.ensembl.org")
}

# See https://www.ensembl.org/info/data/biomart/biomart_restful.html for more details.
get_biomart_xml <- function(dataset) {
stringr::str_glue("
<?xml version='1.0' encoding='UTF-8'?>
<!DOCTYPE Query>
<Query virtualSchemaName = 'default' uniqueRows = '1' count='0' datasetConfigVersion='0.6' header='1' formatter='TSV' requestid='biomaRt'>
  <Dataset name = '{dataset}'>
    <Attribute name = 'ensembl_gene_id'/>
    <Attribute name = 'external_gene_name'/>
    <Attribute name = 'go_id'/>
    <Attribute name = 'go_linkage_type'/>
  </Dataset>
</Query>") |>
    stringr::str_replace_all("\n", "") |>
    stringr::str_replace_all("\\>\\s+\\<", "><")
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
  key <- term_id <- value <- name <- term_name <- namespace <- NULL

  terms <- parsed |>
    dplyr::filter(key %in% c("name", "namespace")) |>
    tidyr::pivot_wider(id_cols = term_id, names_from = key, values_from = value) |>
    dplyr::rename(term_name = name)

  alt_terms <- parsed |>
    dplyr::filter(key == "alt_id") |>
    dplyr::left_join(terms, by = "term_id") |>
    dplyr::select(term_id = value, term_name, namespace)

  dplyr::bind_rows(
    terms,
    alt_terms
  )
}


#' Download GO term descriptions
#'
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with term_id and term_name.
#' @noRd
fetch_go_terms <- function(use_cache, on_error) {
  if(!assert_url_path(get_go_obo_dir(), on_error))
    return(NULL)

  obo_file <- get_go_obo_file()
  lpath <- cached_url_path("obo", obo_file, use_cache)
  readr::read_lines(lpath) |>
    parse_obo_file() |>
    extract_obo_terms()
}

#' Parse GO species YAML
#'
#' @param goex GO species YAML content as a character vector.
#'
#' @return A tibble with GO species metadata.
#' @noRd
parse_go_species_yaml <- function(goex) {
  starts <- stringr::str_which(goex, "^\\s*-\\s+taxon_id:")
  ends <- c(starts[-1] - 1, length(goex))

  purrr::map2(starts, ends, function(i1, i2) {
    organism <- goex[i1:i2]
    organism <- organism[stringr::str_detect(organism, "^\\s*(-\\s+)?[[:alnum:]_]+:")]
    organism <- stringr::str_remove(organism, "^\\s*-\\s+")
    organism <- stringr::str_trim(organism)

    entries <- stringr::str_split_fixed(organism, ":\\s*", 2)
    values <- stringr::str_trim(entries[, 2])
    values[values == ""] <- NA_character_
    stats::setNames(as.list(values), entries[, 1])
  }) |>
    purrr::map(tibble::as_tibble) |>
    purrr::list_rbind()
}

#' Find all species available from geneontology.org
#'
#' This function downloads and parses Gene Ontology's organism metadata table
#' to return the current annotation file species designations. By default, one
#' preferred designation is returned for each species. Model Organism Database
#' (MOD) annotations are preferred when they are available; otherwise UniProt
#' annotations are returned. Set \code{all_sources = TRUE} to return both MOD
#' and UniProt designations.
#'
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#' @param all_sources Logical, if TRUE, return all available annotation sources
#'   for each species. If FALSE, return one preferred designation per species.
#'
#' @return A tibble with Gene Ontology species metadata. Column
#'   \code{designation} contains species designations used in function
#'   \code{fetch_go}.
#' @export
#' @examples
#' go_species <- fetch_go_species(on_error = "warn")
#' go_species_all <- fetch_go_species(on_error = "warn", all_sources = TRUE)
fetch_go_species <- function(on_error = c("stop", "warn", "ignore"),
                             all_sources = FALSE) {
  on_error <- match.arg(on_error)
  assertthat::assert_that(assertthat::is.flag(all_sources))

  # Binding variables from non-standard evaluation locally
  taxon_id <- taxonomic_group <- full_name <- common_name_goc <- code_uniprot <- NULL
  group <- designation <- tax_id <- source_priority <- species <- common_name <- NULL
  annotation_source <- NULL

  resp <- http_request(get_go_species_url(), "")
  if(resp$is_error)
    return(catch_error("Gene Ontology", resp, on_error))

  goex <- httr2::resp_body_string(resp$response) |>
    I() |>
    readr::read_lines()

  go_species <- parse_go_species_yaml(goex) |>
    dplyr::select(taxon_id, taxonomic_group, species = full_name,
                  common_name = common_name_goc, code_uniprot, group) |>
    dplyr::filter(!is.na(code_uniprot))

  go_species_uniprot <- go_species |>
    dplyr::mutate(designation = stringr::str_glue("{code_uniprot}-uniprot"),
                  annotation_source = "UniProt",
                  source_priority = 2)

  go_species_mod <- go_species |>
    dplyr::filter(group != "UniProt") |>
    dplyr::mutate(designation = stringr::str_glue("{code_uniprot}-mod"),
                  annotation_source = "MOD",
                  source_priority = 1)

  go_species <- dplyr::bind_rows(go_species_mod, go_species_uniprot) |>
    dplyr::mutate(tax_id = stringr::str_remove(taxon_id, "^NCBITaxon:")) |>
    dplyr::select(designation, species, common_name, tax_id, taxonomic_group,
                  group, annotation_source, source_priority) |>
    dplyr::arrange(designation) |>
    dplyr::distinct()

  if(!all_sources) {
    go_species <- go_species |>
      dplyr::arrange(tax_id, source_priority, designation) |>
      dplyr::group_by(tax_id) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup()
  }

  go_species |>
    dplyr::select(-source_priority)
}


#' Check if GO species are valid
#'
#' Checks a species argument against current Gene Ontology designations and
#' packaged legacy names. Error messages refer users to exported helper
#' functions.
#'
#' @param species A string, species designation or legacy Gene Ontology species
#'   name.
#' @param on_error A character string indicating the error handling strategy.
#'
#' @return A tibble with valid current and legacy species designations, or
#'   \code{NULL} if current designations cannot be downloaded and \code{on_error}
#'   is not \code{"stop"}.
#' @noRd
assert_go_species <- function(species, on_error) {
  assert_that(is.string(species))

  designation <- legacy_name <- NULL

  go_species <- fetch_go_species(on_error, all_sources = TRUE)
  if(is.null(go_species))
    return(NULL)

  valid_species <- go_species |>
    dplyr::select(designation) |>
    dplyr::bind_rows(go_legacy_mapping |> dplyr::distinct(designation = legacy_name))

  assert_that(species %in% valid_species$designation,
    msg = stringr::str_glue("Invalid species {species}. Use fetch_go_species() to find current species designations or get_go_legacy_mapping() to inspect legacy names.")
  )
  valid_species
}


#' Get Gene Ontology legacy species mapping
#'
#' Returns the packaged mapping between legacy Gene Ontology annotation species
#' names and current annotation file designations. This table is used
#' internally to keep legacy \code{species} values working in \code{fetch_go()},
#' but is also exposed so users can inspect or update their own code. Some
#' legacy names map to more than one current designation.
#'
#' @return A tibble with columns \code{legacy_name} and \code{designation}.
#' @export
#' @examples
#' go_legacy_mapping <- get_go_legacy_mapping()
get_go_legacy_mapping <- function() {
  go_legacy_mapping
}


#' Download GO term gene mapping from geneontology.org
#'
#' @param species Species designation. Base file name for species file under
#'   \url{http://current.geneontology.org/annotations/gaf}. Examples are
#'   \file{HUMAN-uniprot} for human, \file{MOUSE-mod} for mouse, or
#'   \file{YEAST-mod} for yeast.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with columns \code{gene_symbol}, \code{gene_id},
#'   \code{db_id}, \code{term_id}, and \code{evidence}.
#' @noRd
fetch_go_genes_go <- function(species, use_cache, on_error) {
  # Binding variables from non-standard evaluation locally
  gene_id <- db_object_synonym <- symbol <- NULL
  db_id <- go_term <- evidence <- NULL

  url <- get_go_annotation_url()
  gaf_file <- stringr::str_glue("{url}/gaf/{species}.gaf.gz")
  if(!assert_url_path(gaf_file, on_error))
    return(NULL)

  lpath <- cached_url_path(stringr::str_glue("go_gaf_{species}"), gaf_file, use_cache)
  readr::read_tsv(lpath, comment = "!", quote = "", col_names = GAF_COLUMNS,
                  col_types = GAF_TYPES) |>
    dplyr::mutate(gene_id = stringr::str_remove(db_object_synonym, "\\|.*$")) |>
    dplyr::select(gene_symbol = symbol, gene_id, db_id, term_id = go_term, evidence) |>
    dplyr::distinct()
}


#' Match legacy GO species names to current designations
#'
#' @param species Species designation or legacy Gene Ontology species name.
#' @param on_error A character string indicating the error handling strategy.
#'
#' @return A current Gene Ontology species designation, or \code{NULL} if
#'   \code{species} is an ambiguous legacy name and \code{on_error} is not
#'   \code{"stop"}.
#' @noRd
match_go_legacy_species <- function(species, on_error) {
 if(species %in% go_legacy_mapping$legacy_name) {
    new_species <- go_legacy_mapping$designation[go_legacy_mapping$legacy_name == species]
    if(length(new_species) > 1) {
      des <- stringr::str_c(new_species, collapse = ', ')
      msg <- stringr::str_glue("Legacy species {species} corresponds to multiple current designations. Please choose one of: {des}")
      return(error_response(msg, on_error))
    }
    species <- new_species
  }
  species
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
#'   Object Synonym) is returned as \code{gene_id}. It is up to the user to
#'   select the appropriate database.
#'
#' @param species Species designation. Examples are \file{HUMAN-uniprot} for
#'   human, \file{MOUSE-mod} for mouse, or \file{YEAST-mod} for yeast. Legacy
#'   names such as \file{goa_human}, \file{mgi}, and \file{sgd} are also
#'   accepted when they map to a single current designation. Current species
#'   designations can be obtained using \code{fetch_go_species}; legacy mappings
#'   can be inspected using \code{get_go_legacy_mapping}.
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @importFrom assertthat assert_that
#' @noRd
fetch_go_from_go <- function(species, use_cache, on_error) {
  assert_that(!missing(species), msg = "Argument 'species' is missing.")
  assert_go_species(species, on_error)

  species <- match_go_legacy_species(species, on_error)
  if(is.null(species))
    return(NULL)
 
  mapping <- fetch_go_genes_go(species = species, use_cache = use_cache, on_error = on_error)
  if(is.null(mapping))
    return(NULL)

  terms <- fetch_go_terms(use_cache = use_cache, on_error = on_error)
  if(is.null(terms))
    return(NULL)

  list(
    terms = terms,
    mapping = mapping
  )
}



#' Download GO term gene mapping from Ensembl
#'
#' @param dataset Dataset you want to use. To see the different datasets
#'   available within a biomaRt you can e.g. do: mart <-
#'   biomaRt::useEnsembl(biomart = "ensembl"), followed by
#'   biomaRt::listDatasets(mart).
#' @param use_cache Logical, if TRUE, the remote data will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A tibble with columns \code{gene_id}, \code{gene_symbol},
#'   \code{term_id} and \code{evidence}.
#' @noRd
fetch_go_genes_bm <- function(dataset, use_cache, on_error) {
  xml <- get_biomart_xml(dataset) |>
    stringr::str_replace_all("\\s", "%20")
  biomart_path <- paste0(get_biomart_url(), "/biomart/martservice")
  if(!assert_url_path(biomart_path, on_error))
    return(NULL)

  req <- paste0(biomart_path, "?query=", xml)

  # Problems with cache, bfcneedsupdate returns error for this query
  # lpath <- cached_url_path(stringr::str_glue("biomart_{dataset}"), resp, use_cache)
  res <- readr::read_tsv(req, show_col_types = FALSE)
  if(ncol(res) == 4) {
    res |> rlang::set_names(c("gene_id", "gene_symbol", "term_id", "evidence"))
  } else {
    error_response("Problem with Biomart", on_error)
  }
}



#' Get functional term data from Ensembl
#'
#' Download term information (GO term ID and name) and gene-term mapping
#' (gene ID, symbol and GO term ID) from Ensembl.
#'
#' @param dataset Dataset you want to use. To see the different datasets
#'   available within a biomaRt you can e.g. do: mart <-
#'   biomaRt::useEnsembl(biomart = "ensembl"), followed by
#'   biomaRt::listDatasets(mart).
#' @param use_cache Logical, if TRUE, the remote file will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @importFrom assertthat assert_that is.string
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @noRd
fetch_go_from_bm <- function(dataset, use_cache, on_error) {
  assert_that(!missing(dataset), msg = "Argument 'dataset' is missing.")
  assert_that(is.string(dataset))

  mapping <- fetch_go_genes_bm(dataset, use_cache = use_cache, on_error = on_error)
  if(is.null(mapping))
    return(error_response("Could not retrieve mapping from Ensembl", on_error))
  terms <- fetch_go_terms(use_cache = use_cache, on_error)

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
#'   Synonym) is returned as \code{gene_id}. It is up to the user to select
#'   the appropriate database.
#'
#'   Alternatively, if \code{dataset} is provided, mapping will be downloaded
#'   from Ensembl database. It will contain gene symbols and Ensembl gene IDs.
#'
#' @param species (Optional) Species designation. Examples are
#'   \code{HUMAN-uniprot} for human, \code{MOUSE-mod} for mouse, or
#'   \code{YEAST-mod} for yeast. Legacy names such as \code{goa_human},
#'   \code{mgi}, and \code{sgd} are also accepted when they map to a single
#'   current designation. Current species designations can be obtained using
#'   \code{fetch_go_species}; legacy mappings can be inspected using
#'   \code{get_go_legacy_mapping}. This argument is used when fetching data from
#'   the Gene Ontology database.
#' @param dataset (Optional) A string representing the dataset passed to
#'   Ensembl's Biomart, e.g. 'scerevisiae_gene_ensembl'. To see the different
#'   datasets available within a biomaRt you can e.g. do: mart <-
#'   biomaRt::useEnsembl(biomart = "ensembl"), followed by
#'   biomaRt::listDatasets(mart).
#' @param use_cache Logical, if TRUE, the remote data will be cached locally.
#' @param on_error A character string indicating the error handling strategy:
#'   either "stop" to halt execution, "warn" to issue a warning and return
#'   `NULL` or "ignore" to return `NULL` without warnings. Defaults to "stop".
#'
#' @return A list with \code{terms} and \code{mapping} tibbles.
#' @export
#' @importFrom assertthat assert_that
#' @examples
#' # Fetch GO data from Ensembl
#' go_data_ensembl <- fetch_go(dataset = "scerevisiae_gene_ensembl", on_error = "warn")
#' # Fetch GO data from Gene Ontology
#' go_data_go <- fetch_go(species = "YEAST-mod", on_error = "warn")
fetch_go <- function(species = NULL, dataset = NULL, use_cache = TRUE,
                     on_error = c("stop", "warn", "ignore")) {
  on_error <- match.arg(on_error)

  assert_that(!(is.null(species) & is.null(dataset)),
              msg = "One of the arguments 'species' or 'dataset' must be supplied.")
  assert_that(is.null(species) | is.null(dataset),
              msg = "Only one of the arguments 'species' or 'dataset' must be supplied.")

  if (!is.null(species)) {
    fetch_go_from_go(species, use_cache = use_cache, on_error = on_error)
  } else {
    fetch_go_from_bm(dataset, use_cache = use_cache, on_error = on_error)
  }
}
