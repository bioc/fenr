library(testthat)


# Set 100 features
N <- 100
features_all <- sprintf("gene_%03d", seq_len(N))

# Three terms with 50, 10 and 3 features
term_ids <- c("term_1", "term_2", "term_3")
term_names <- c("name_1", "name_2", "name_3")
term_sizes <- c(50, 10, 3)


terms <- tibble::tibble(
  term_id = term_ids,
  term_name = term_names
)

# random selection of features for terms
set.seed(666)
mapping <- purrr::map2_dfr(term_ids, term_sizes, function(tid, n) {
  tibble::tibble(
    term_id = tid,
    feature_id = sample(features_all, n)
  )
}) |>
  dplyr::add_row(term_id = "missing_term", feature_id = "gene_000")

all_terms <- mapping$term_id |> unique()
all_features <- mapping$feature_id |> unique()


test_that("Expected correct output", {
  td <- prepare_for_enrichment(terms, mapping, feature_name = "feature_id")

  # Correct output
  expect_equal(length(td), 3)

  # Check class
  expect_s3_class(td, "fenr_terms")

  # Check term names
  expect_equal(length(td$term2name), length(all_terms))
  for(i in seq_along(terms$term_id)) {
    r <- terms[i, ]
    expect_equal(sort(r$term_name), sort(td$term2name[[r$term_id]]))
  }

  # Check for a missing term
  expect_contains(names(td$term2name), "missing_term")

  # Check term-feature hash
  term_ids <- mapping$term_id |> unique()
  chk1 <- term_ids |>
    purrr::map(function(trm) {
      expected <- mapping |>
        dplyr::filter(term_id == trm) |>
        dplyr::pull(feature_id) |>
        sort()
      returned <- td$term2feature[[trm]] |>
        sort()
      returned_fun <- get_term_features(td, trm) |>
        sort()
      expect_equal(expected, returned)
      expect_equal(expected, returned_fun)
    })

  # Check feature-term hash
  expect_equal(length(td$feature2term), length(all_features))
  chk2 <- all_features |>
    purrr::map(function(feat) {
      expected <- mapping |>
        dplyr::filter(feature_id == feat) |>
        dplyr::pull(term_id) |>
        sort()
      returned <- td$feature2term[[feat]] |>
        sort()
      returned_fun <- get_feature_terms(td, feat) |>
        sort()
      expect_equal(expected, returned)
      expect_equal(expected, returned_fun)
    })
})


test_that("Duplicate term-gene pairs are removed", {
  extra_row <- mapping[1, ]
  mapping_dup <- mapping |>
    dplyr::add_row(extra_row)
  td <- prepare_for_enrichment(terms, mapping_dup, feature_name = "feature_id")
  returned <- td$term2feature[[extra_row$term_id]] |>
    sort()
  expected <- mapping |>
    dplyr::filter(term_id == extra_row$term_id) |>
    dplyr::pull(feature_id) |>
    sort()
  expect_equal(expected, returned)
})


test_that("Wrong columns in terms", {
  trms <- terms |>
    dplyr::rename(gobble = term_id)

  expect_error(
    prepare_for_enrichment(trms, mapping, feature_name = "feature_id")
  )
})


test_that("Missing column in mapping", {
  mp <- mapping |>
    dplyr::rename(mobble = term_id)

  expect_error(
    prepare_for_enrichment(terms, mp, feature_name = "feature_id")
  )
})


test_that("Wrong feature name", {
  expect_error(
    prepare_for_enrichment(terms, mapping, feature_name = "trouble")
  )
})


test_that("No overlap between all_features and mapping", {
  expect_error(
    prepare_for_enrichment(terms, mapping, all_features = c("a", "b", "c"),  feature_name = "feature_id")
  )
})
