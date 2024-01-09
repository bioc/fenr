library(testthat)

# Set 100 features
N <- 100
features_all <- sprintf("gene_%03d", seq_len(N))
symbols_all <- sprintf("symbol_%03d", seq_len(N))
feat2name <- purrr::set_names(symbols_all, features_all)

# Three terms with 50, 10 and 3 features
term_ids <- c("term_1", "term_2", "term_3")
term_names <- c("name_1", "name_2", "name_3")
term_sizes <- c(50, 10, 3)

# Prepare data for enrichment. The last term's name is intentionally missing.
term2name <- new.env(hash = TRUE)
for (i in seq_len(length(term_ids) - 1)) {
  term2name[[term_ids[i]]] <- term_names[i]
}

# random selection of features for terms
set.seed(666)
mapping <- purrr::map2_dfr(term_ids, term_sizes, function(tid, n) {
  tibble::tibble(
    term_id = tid,
    feature_id = sample(features_all, n)
  )
})

# Feature to terms hash
feature2term <- new.env(hash = TRUE)
features <- unique(mapping$feature_id)
for(feat in features) {
  terms <- mapping |>
    dplyr::filter(feature_id == feat) |>
    dplyr::pull(term_id)
  feature2term[[feat]] <- terms
}

# Term to feature hash
term2feature <- new.env(hash = TRUE)
terms <- unique(mapping$term_id)
for(term in terms) {
  features <- mapping |>
    dplyr::filter(term_id == term) |>
    dplyr::pull(feature_id)
  term2feature[[term]] <- features
}

# final structure required by functional_enrichment
term_data <- list(
  term2name = term2name,
  feature2term = feature2term,
  term2feature = term2feature
) |>
  structure(class = "fenr_terms")

term_stats <- function(tid, N_sel, n_with_sel) {
  N_with <- length(term2feature[[tid]])
  n_without_sel <- N_sel - n_with_sel
  n_expect <- N_with * N_sel / N
  enrichment <- n_with_sel / n_expect

  n_with_nsel <- N_with - n_with_sel
  n_without_nsel <- N - (n_with_sel + n_without_sel + n_with_nsel)
  cont_tab <- rbind(
    c(n_with_sel, n_with_nsel),
    c(n_without_sel, n_without_nsel)
  )

  # in the test we use fisher.test function as opposed to phyper in the code
  # p <- 1 - stats::phyper(n_with_sel - 1, N_with, N_without, N_sel)
  p <- fisher.test(cont_tab, alternative = "greater")$p.value

  sel_with <- term2feature[[tid]][seq_len(n_with_sel)]
  if (n_without_sel > 0) {
    sel_without <- setdiff(features_all, term2feature[[tid]])[seq_len(n_without_sel)]
  } else {
    sel_without <- character(0)
  }
  sel <- c(sel_with, sel_without)

  list(
    N_with = N_with,
    n_with_sel = n_with_sel,
    n_expect = n_expect,
    enrichment = enrichment |> signif(3),
    p = p |> signif(3),
    sel = sel
  )
}

# Test selection of 3 genes
set.seed(42)
test_sel <- mapping |>
  dplyr::group_by(term_id) |>
  dplyr::sample_n(3) |>
  dplyr::pull(feature_id)


#######################################################

test_that("Incorrect term_data class", {
  td <- term_data
  class(td) <- "wobble"

  expect_error(
    functional_enrichment(features_all, features_all[seq_len(10)], td)
  )
})


test_that("Features must be character vectors with length > 1", {
  expect_error(functional_enrichment(c(1, 2, 3, 4), c(1, 2), term_data))
  expect_error(functional_enrichment(features_all, features_all[1], term_data))
  expect_error(functional_enrichment(features_all, character(0), term_data))
  expect_error(functional_enrichment(features_all[1], features_all[1], term_data))
  expect_error(functional_enrichment(character(0), character(0), term_data))
})

test_that("Correct output", {
  enr <- functional_enrichment(features_all, test_sel, term_data)

  expect_true(is.data.frame(enr))
  expect_equal(nrow(enr), length(term_ids))
  expect_equal(ncol(enr), 10)
})

test_that("Gene names translated", {
  enr <- functional_enrichment(features_all, test_sel, term_data, feat2name)
  symbols <- enr |>
    tidyr::separate_longer_delim(ids, delim = ", ") |>
    dplyr::pull(ids) |>
    unique()
  expect_contains(symbols_all, symbols)
})


test_that("Missing term name replaced with NA", {
  enr <- functional_enrichment(features_all, test_sel, term_data)
  misname <- enr |>
    dplyr::filter(term_id == "term_3") |>
    dplyr::pull(term_name)
  expect_true(is.na(misname))
})


test_that("n_with_sel < 2 should return NULL", {
  set.seed(42)
  # genes for term 1
  feat_1 <- mapping |>
    dplyr::filter(term_id == term_ids[1]) |>
    dplyr::pull(feature_id)
  # genes for term 2
  feat_2 <- mapping |>
    dplyr::filter(term_id == term_ids[2]) |>
    dplyr::pull(feature_id)
  # find genes exclusive to selections
  sel_1 <- setdiff(feat_1, feat_2)
  sel_2 <- setdiff(feat_2, feat_1)

  # One gene from feat_1, 3 genes from feat_2
  sel <- c(sample(sel_1, 1), sample(sel_2, 3))
  enr <- functional_enrichment(features_all, sel, term_data)
  # Term 1 should not appear in the output
  expect_false(term_ids[1] %in% enr$term_id)

  # One gene from feat_1, 1 gene from feat_2
  sel <- c(sample(sel_1, 1), sample(sel_2, 1))
  enr <- functional_enrichment(features_all, sel, term_data)
  # All searches have n_with_sel < 2, expect NULL
  expect_null(enr)
})


test_that("Enrichment of 1", {
  tid <- term_ids[1]
  N_sel <- 20
  n_with_sel <- 10

  res <- term_stats(tid, N_sel, n_with_sel)
  enr <- functional_enrichment(features_all, res$sel, term_data) |>
    dplyr::filter(term_id == tid)

  expect_equal(
    c(res$N_with, res$n_with_sel, res$n_expect, res$enrichment, res$p),
    c(enr$N_with, enr$n_with_sel, enr$n_expect, enr$enrichment, enr$p_value)
  )
})


test_that("All selected features in term", {
  tid <- term_ids[2]
  N_sel <- term_sizes[2]
  n_with_sel <- term_sizes[2]

  res <- term_stats(tid, N_sel, n_with_sel)
  enr <- functional_enrichment(features_all, res$sel, term_data) |>
    dplyr::filter(term_id == tid)

  expect_equal(
    c(res$N_with, res$n_with_sel, res$n_expect, res$enrichment, res$p),
    c(enr$N_with, enr$n_with_sel, enr$n_expect, enr$enrichment, enr$p_value)
  )
})


test_that("No selected features in term", {
  tid <- "term_3"
  N_sel <- 20
  n_with_sel <- 0

  res <- term_stats(tid, N_sel, n_with_sel)
  enr <- functional_enrichment(features_all, res$sel, term_data) |>
    dplyr::filter(term_id == tid)

  expect_equal(nrow(enr), 0)
})


test_that("Selection is all features", {
  enr <- functional_enrichment(features_all, features_all, term_data)

  n <- nrow(enr)
  expect_equal(enr$enrichment, rep(1, n))
  expect_equal(enr$p_value, rep(1, n))
})


test_that("No match between selection and all features", {
  expect_null(
    functional_enrichment(features_all, c("a", "b", "c"), term_data)
  )
})


