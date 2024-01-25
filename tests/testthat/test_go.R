test_that("Incorrect aruments in fetch_go", {
  expect_error(fetch_go())
  expect_error(fetch_go(species = "sgd", mart = "mart"))
})


test_that("Incorrect species in fetch_go", {
  expect_error(fetch_go(species = 1243))
  expect_error(fetch_go(species = "not a species"))
})


test_that("Incorrect mart in fetch_go", {
  expect_error(fetch_go(mart = 1))
  expect_error(fetch_go(mart = "mart"))
})


test_that("Processing OBO file", {
  expected_ids <- c("GO:0000001", "GO:0000022", "GO:0000232", "GO:1905121")
  expected_names <- c("mitochondrion inheritance", "mitotic spindle elongation",
                      "obsolete nuclear interphase chromosome", "mitotic spindle elongation")

  trms <- readr::read_lines("../test_data/go_obo_test.txt") |>
    parse_obo_file() |>
    extract_obo_terms()

  expect_s3_class(trms, "tbl")
  expect_equal(trms$term_id, expected_ids)
  expect_equal(trms$term_name, expected_names)
})


test_that("Expected behaviour from a non-responsive server", {
  species <- "sgd"
  httr2::with_mocked_responses(
    mock = mocked_500,
    code = {
      test_unresponsive_server(fetch_go_terms, use_cache = FALSE)
      test_unresponsive_server(fetch_go_species)
      test_unresponsive_server(fetch_go_from_go, species = species, use_cache = FALSE)
    })
})


test_that("Expected return from fetch_go_species", {
  expected_selection <- c("goa_human", "mgi", "rgd", "sgd", "fb", "wb")
  spec <- fetch_go_species(on_error = "warn")
  if(!is.null(spec)) {
    expect_is(spec, "tbl")
    expect_true(all(expected_selection %in% spec$designation))
  }
})


test_that("GO yeast from GO is correct", {
  species <- "sgd"

  expected_terms <- tibble::tribble(
    ~term_id, ~term_name,
    "GO:0000166", "nucleotide binding",
    "GO:0006096", "glycolytic process",
    "GO:0005199", "structural constituent of cell wall",
    "GO:0004365", "glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity"
  )

  expected_mapping <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "GO:0000166", "POL1",
    "GO:0006096", "FBA1",
    "GO:0005199", "CWP2",
    "GO:0004365", "TDH3"
  )

  re <- fetch_go(species = species, on_error = "warn")
  if(!is.null(re)) {
    test_fetched_structure(re)
    test_terms(re$terms, expected_terms)
    test_mapping(re$mapping, expected_mapping, "gene_symbol")
  }
})


test_that("GO yeast from Ensembl is correct", {
  expected_terms <- tibble::tribble(
    ~term_id, ~term_name,
    "GO:0000166", "nucleotide binding",
    "GO:0006096", "glycolytic process",
    "GO:0005199", "structural constituent of cell wall",
    "GO:0004365", "glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity"
  )

  expected_mapping <- tibble::tribble(
    ~term_id, ~ensembl_gene_id,
    "GO:0000166", "YNL102W",
    "GO:0006096", "YKL060C",
    "GO:0005199", "YKL096W-A",
    "GO:0004365", "YGR192C"
  )

  mart <- tryCatch(
    {biomaRt::useEnsembl(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")},
    error = function(cond) {
      message("Biomart is not responding.")
      message(conditionMessage(cond))
      NULL
    }
  )

  if(!is.null(mart)) {
    re <- fetch_go(mart = mart)
    test_fetched_structure(re)
    test_terms(re$terms, expected_terms)
    test_mapping(re$mapping, expected_mapping, "ensembl_gene_id")
  }
})
