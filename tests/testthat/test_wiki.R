test_that("Incorrect species in fetch_wiki", {
  expect_error(fetch_wiki())
  expect_error(fetch_wiki(1243))
  expect_error(fetch_wiki("not a species"))
})


test_that("Incorrect species in fetch_wiki_pathways", {
  expect_error(fetch_wiki_pathways("Heffalump"))
})


test_that("Parsing GPLM", {
  expected <- tibble::tibble(
    TextLabel = c("STR3", "CYS4", "STR2"),
    Type = c("GeneProduct", "GeneProduct", "GeneProduct"),
    Database = c("SGD", "SGD", "SGD"),
    ID = c("S000003152", "S000003387", "S000003891")
  )

  gpml <- readLines("../test_data/wiki_test.gpml") |>
    paste(collapse = "\n")
  res <- parse_wiki_gpml(gpml)

  expect_equal(res, expected)
})

test_that("Expected behaviour from a non-responsive server", {
  species <- "Bacillus subtilis"
  pathway <- "WP1466"
  httr2::with_mocked_responses(
    mock = mocked_500,
    code = {
      test_unresponsive_server(fetch_wiki_species)
      test_unresponsive_server(fetch_wiki_pathways, species = species)
      test_unresponsive_server(fetch_wiki, species = species)
      test_unresponsive_server(fetch_wiki_pathway_genes_api, pathways = pathway)
      test_unresponsive_server(fetch_wiki, species = species)
    })
})

test_that("Expected return from fetch_wiki_species", {
  expected_selection <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus",
                          "Saccharomyces cerevisiae", "Drosophila melanogaster", "Caenorhabditis elegans")
  spec <- fetch_wiki_species(on_error = "warn")
  if(!is.null(spec)) {
    expect_is(spec, "tbl")
    expect_true(all(expected_selection %in% spec$designation))
  }
})


test_that("WikiPathways correct response", {
  species <- "Bacillus subtilis"
  databases <- c("Ensembl", "Entrez Gene", "HGNC", "HGNC Accession number", "Uniprot-TrEMBL")
  types <- c("GeneProduct", "Protein", "Rna", "RNA")

  expected_terms <- tibble::tribble(
    ~term_id, ~term_name,
    "WP1466", "Response regulator aspartate phosphatase interactions",
    "WP1527", "Stress response",
    "WP2360", "Folate biosynthesis"
  )

  expected_mapping <- tibble::tribble(
    ~term_id, ~text_label,
    "WP1466", "rapA",
    "WP1466", "rapE",
    "WP1527", "sinR",
    "WP1527", "sinI"
  )

  re <- fetch_wiki(species, databases = databases, types = types, on_error = "warn")
  if(!is.null(re)) {
    test_fetched_structure(re)
    test_terms(re$terms, expected_terms)
    test_mapping(re$mapping, expected_mapping, "text_label")
  }
})
