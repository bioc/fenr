test_that("Incorrect species in fetch_kegg", {
  expect_error(fetch_kegg())
  expect_error(fetch_kegg(1243))
  expect_error(fetch_kegg("not a species"))
})

test_that("Incorrect batch size in fetch_kegg", {
  expect_error(fetch_kegg("sce", "aha!"))
  expect_error(fetch_kegg("sce", -1))
  expect_error(fetch_kegg("sce", 0))
  expect_error(fetch_kegg("sce", 1000))
})


test_that("Parsing KEGG flat file", {
  flat <- readr::read_file("../test_data/kegg_test.txt")
  expected <- readr::read_rds("../test_data/kegg_parsing_result.rds")
  parsed <- parse_kegg_genes(flat)
  expect_equal(parsed, expected)
})


test_that("Expected behaviour from a non-responsive server", {
  species <- "mge"
  pathway <- "mge00010"
  httr2::with_mocked_responses(
    mock = mocked_500,
    code = {
      test_unresponsive_server(fetch_kegg_species)
      test_unresponsive_server(fetch_kegg_pathways, species = species)
      test_unresponsive_server(fetch_kegg_mapping, pathways = pathway, batch_size = 1)
      test_unresponsive_server(fetch_kegg, species = species)
    })
})


test_that("Expected return from fetch_kegg_species", {
  expected_selection <- c("hsa", "mmu", "rno", "sce", "dme", "cel")
  spec <- fetch_kegg_species(on_error = "warn")
  if(!is.null(spec)) {
    expect_is(spec, "tbl")
    expect_true(all(expected_selection %in% spec$designation))
  }
})


test_that("Correct response from fetch_kegg", {
  species <- "mge"

  expected_terms <- tibble::tribble(
    ~term_id, ~term_name,
    "mge01100", "Metabolic pathways - Mycoplasmoides genitalium G37",
    "mge01230", "Biosynthesis of amino acids - Mycoplasmoides genitalium G37",
    "mge03010", "Ribosome - Mycoplasmoides genitalium G37",
    "mge03030", "DNA replication - Mycoplasmoides genitalium G37"
  )

  expected_mapping <- tibble::tribble(
    ~term_id, ~gene_symbol,
    "mge00010", "pgi",
    "mge00010", "fba",
    "mge03010", "rpsB",
    "mge03030", "polC"
  )

  re <- fetch_kegg(species, on_error = "warn")
  if(!is.null(re)) {
    test_fetched_structure(re)
    test_terms(re$terms, expected_terms)
    test_mapping(re$mapping, expected_mapping, "gene_symbol")
  }
})
