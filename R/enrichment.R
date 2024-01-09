#' Prepare Term Data for Enrichment Analysis
#'
#' Process term data downloaded with the \code{fetch_*} functions, preparing it
#' for fast enrichment analysis using \code{functional_enrichment}.
#'
#' @details
#'
#' This function takes two tibbles containing functional term information
#' (\code{terms}) and feature mapping (\code{mapping}), and converts them into
#' an object required by \code{functional_enrichment} for efficient analysis.
#' Terms and mapping can be generated with the database access functions
#' included in this package, such as \code{fetch_reactome} or
#' \code{fetch_go_from_go}.
#'
#' @param terms A tibble with at least two columns: \code{term_id} and
#'   \code{term_name}. This tibble contains information about functional term
#'   names and descriptions.
#' @param mapping A tibble with at least two columns, containing the mapping
#'   between functional terms and features. One column must be named
#'   \code{term_id}, while the other column should have a name specified by the
#'   \code{feature_name} argument. For example, if \code{mapping} contains
#'   columns \code{term_id}, \code{accession_number}, and \code{gene_symbol},
#'   setting \code{feature_name = "gene_symbol"} indicates that gene symbols
#'   will be used for enrichment analysis.
#' @param all_features A vector with all feature IDs used as the background for
#'   enrichment. If not specified, all features found in \code{mapping} will be
#'   used, resulting in a larger object size.
#' @param feature_name The name of the column in the \code{mapping} tibble to be
#'   used as the feature identifier. For example, if \code{mapping} contains
#'   columns \code{term_id}, \code{accession_number}, and \code{gene_symbol},
#'   setting \code{feature_name = "gene_symbol"} indicates that gene symbols
#'   will be used for enrichment analysis.
#'
#' @return An object of class \code{fenr_terms} required by
#'   \code{functional_enrichment}.
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' data(exmpl_all)
#' go <- fetch_go(species = "sgd")
#' go_terms <- prepare_for_enrichment(go$terms, go$mapping, exmpl_all,
#'                                    feature_name = "gene_symbol")
prepare_for_enrichment <- function(terms, mapping, all_features = NULL,
                                   feature_name = "gene_id") {
  # Binding variables from non-standard evaluation locally
  feature_id <- term_id <- NULL

  # Argument checks
  assert_that(is.data.frame(terms) || tibble::is_tibble(terms),
              msg = "'terms' must be a data frame or tibble.")

  assert_that(is.data.frame(mapping) || tibble::is_tibble(mapping),
              msg = "'mapping' must be a data frame or tibble.")

  assert_that(is.null(all_features) || is.vector(all_features),
              msg = "'all_features' must be a vector or NULL.")

  assert_that(is.character(feature_name) && length(feature_name) == 1,
              msg = "'feature_name' must be a single string.")

  # Check terms
  assert_that(all(c("term_id", "term_name") %in% colnames(terms)),
              msg = "Column names in 'terms' should be 'term_id' and 'term_name'.")

  assert_that(anyDuplicated(terms$term_id) == 0,
              msg = "Duplicated term_id detected in 'terms'.")

  # Check mapping
  assert_that("term_id" %in% colnames(mapping),
              msg = "'mapping' should contain a column named 'term_id'.")

  # Check for feature name
  assert_that(feature_name %in% colnames(mapping),
              msg = paste0(feature_name, " column not found in mapping table.
                           Check 'feature_name' argument."))

  # Replace empty all_features with everything from mapping
  map_features <- mapping[[feature_name]] |>
    unique()
  if (is.null(all_features)) {
    all_features <- map_features
  } else {
    # Check if mapping is contained in all features
    if (length(intersect(all_features, map_features)) == 0)
      stop("No overlap between 'all_features' and features found in 'mapping'.
           Did you provide correct 'all_features'?")
  }

  # Check for missing term descriptions
  mis_term <- setdiff(mapping$term_id, terms$term_id)
  if (length(mis_term) > 0) {
    dummy <- tibble::tibble(
      term_id = mis_term,
      term_name = rep(NA_character_, length(mis_term))
    )
    terms <- dplyr::bind_rows(terms, dummy)
  }

  # Hash to select term name
  term2name <- new.env(hash = TRUE)
  for (i in seq_len(nrow(terms))) {
    r <- terms[i, ]
    term2name[[r$term_id]] <- r$term_name
  }

  # feature-term tibble
  feature_term <- mapping |>
    dplyr::rename(feature_id = !!feature_name) |>
    dplyr::filter(feature_id %in% all_features) |>
    dplyr::select(feature_id, term_id)

  # Feature to terms hash
  f2t <- feature_term |>
    dplyr::group_by(feature_id) |>
    dplyr::summarise(terms = list(term_id)) |>
    tibble::deframe()
  feature2term <- new.env(hash = TRUE)
  for(feat in names(f2t))
    feature2term[[feat]] <- f2t[[feat]]

  # Term to feature hash
  t2f <- feature_term |>
    dplyr::group_by(term_id) |>
    dplyr::summarise(features = list(feature_id)) |>
    tibble::deframe()
  term2feature <- new.env(hash = TRUE)
  for(term in names(t2f))
    term2feature[[term]] <- t2f[[term]]

  list(
    term2name = term2name,
    term2feature = term2feature,
    feature2term = feature2term
  ) |>
    structure(class = "fenr_terms")
}



#' Fast Functional Enrichment
#'
#' Perform fast functional enrichment analysis based on the hypergeometric
#' distribution. Designed for use in interactive applications.
#'
#' @details This function carries out functional enrichment analysis on a
#'   selection of features (e.g., differentially expressed genes) using the
#'   hypergeometric probability distribution (Fisher's exact test). Features can
#'   be genes, proteins, etc. The \code{term_data} object contains functional
#'   term information and feature-term mapping.
#'
#' @param feat_all A character vector with all feature identifiers, serving as
#'   the background for enrichment.
#' @param feat_sel A character vector with feature identifiers in the selection.
#' @param term_data An object of class \code{fenr_terms}, created by
#'   \code{prepare_for_enrichment}.
#' @param feat2name An optional named list to convert feature IDs into feature
#'   names.
#'
#' @return A tibble with enrichment results, providing the following information
#'   for each term:
#'   \itemize{
#'     \item{\code{N_with} - number of features with this term among all features}
#'     \item{\code{n_with_sel} - number of features with this term in the selection}
#'     \item{\code{n_expect} - expected number of features with this term in the selection,
#'       under the null hypothesis that terms are mapped to features randomly}
#'     \item{\code{enrichment} - ratio of n_with_sel / n_expect}
#'     \item{\code{odds_ratio} - odds ratio for enrichment; is infinite when all
#'       features with the given term are in the selection}
#'     \item{\code{p_value} - p-value from a single hypergeometric test}
#'     \item{\code{p_adjust} - p-value adjusted for multiple tests using the
#'       Benjamini-Hochberg approach}
#'   }.
#'
#' @importFrom assertthat assert_that
#' @importFrom methods is
#' @export
#' @examples
#' data(exmpl_all, exmpl_sel)
#' go <- fetch_go(species = "sgd")
#' go_terms <- prepare_for_enrichment(go$terms, go$mapping, exmpl_all, feature_name = "gene_symbol")
#' enr <- functional_enrichment(exmpl_all, exmpl_sel, go_terms)
functional_enrichment <- function(feat_all, feat_sel, term_data, feat2name = NULL) {

  # Binding variables from non-standard evaluation locally
  N_with <- n_with_sel <- n_expect <- enrichment <- odds_ratio <- NULL
  desc <- p_value <- p_adjust <- NULL

  # Check for character vectors
  assert_that(is.character(feat_all))
  assert_that(is.character(feat_sel))
  assert_that(length(feat_all) > 1)
  assert_that(length(feat_sel) > 1)

  # Check term_data class
  assert_that(is(term_data, "fenr_terms"))

  # If no overlap between selection and all, return NULL
  if(!any(feat_sel %in% feat_all))
    return(NULL)

  # all terms present in the selection
  our_terms <- feat_sel |>
    purrr::map(~term_data$feature2term[[.x]]) |>
    unlist() |>
    unique()

  # number of features in selection
  N_sel <- length(feat_sel)
  # total number of features
  N_tot <- length(feat_all)

  res <- purrr::map_dfr(our_terms, function(term_id) {
    # all features with the term
    # term_data$term2feature is a hash environment
    tfeats <- term_data$term2feature[[term_id]]
    # necessary if term data contain features not present in feat_all
    tfeats <- tfeats[tfeats %in% feat_all]

    # features from selection with the term
    # this is faster than intersect(tfeats, feat_sel)
    tfeats_sel <- tfeats[tfeats %in% feat_sel]

    N_with <- length(tfeats)
    N_without <- N_tot - N_with

    # building contingency table
    n_with_sel <- length(tfeats_sel)
    n_without_sel <- N_sel - n_with_sel
    n_with_nsel <- N_with - n_with_sel
    n_without_nsel <- N_tot - (n_with_sel + n_without_sel + n_with_nsel)

    # contingency table is
    #               | In selection  | Not in selection
    #-------------------------------------------------
    #  With term    | n_with_sel    | n_with_nsel
    #  Without term | n_without_sel | n_without_nsel

    if (n_with_sel < 2) return(NULL)

    # Expected number of features in selection, if random
    n_expect <- N_with * N_sel / N_tot

    # Odds ratio
    odds_ratio <- (n_with_sel / n_without_sel) / (n_with_nsel / n_without_nsel)

    # Hypergeometric function much faster than fisher.test
    p <- 1 - stats::phyper(n_with_sel - 1, N_with, N_without, N_sel)

    # Convert feature IDs to feature names;
    if (!is.null(feat2name))
      tfeats_sel <- feat2name[tfeats_sel] |> unname()

    term_name <- term_data$term2name[[term_id]]
    # returns NAs if no term found
    if (is.null(term_name)) term_name <- NA_character_

    # constructing a tibble in every iteration is more time expensive than `c`,
    # even with overhead of converting types afterwards. Not elegant, but fast.
    c(
      term_id = term_id,
      term_name = term_name,
      N_with = N_with,
      n_with_sel = n_with_sel,
      n_expect = n_expect,
      enrichment = n_with_sel / n_expect,
      odds_ratio = odds_ratio,
      ids = paste(tfeats_sel, collapse = ", "),
      p_value = p
    )
  })

  # Drawback - if all selections below minimum, res is tibble 0 x 0, need to
  # catch it
  if (nrow(res) == 0) {
    res <- NULL
  } else {
    res <- res |>
      dplyr::mutate(
        dplyr::across(c(N_with, n_with_sel), as.integer),
        dplyr::across(c(n_expect, enrichment, odds_ratio, p_value), as.numeric),
        p_adjust = stats::p.adjust(p_value, method = "BH"),
        dplyr::across(c(enrichment, odds_ratio, p_value, p_adjust), ~signif(.x, 3)),
        n_expect = round(n_expect, 2)
      ) |>
      dplyr::arrange(desc(odds_ratio))
  }
  res
}
