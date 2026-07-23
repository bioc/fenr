test_that("Incorrect aruments in fetch_go", {
  expect_error(fetch_go())
  expect_error(fetch_go(species = "YEAST-mod", dataset = "scerevisiae_gene_ensembl"))
})


test_that("Incorrect species in fetch_go", {
  expect_error(fetch_go(species = 1243))
  expect_error(fetch_go(species = "not a species"))
})

test_that("Incorrect dataset triggers error", {
  expect_error(fetch_go(dataset = "not a dataset"))
  expect_error(fetch_go_from_bm(dataset = "not a dataset", use_cache = FALSE, on_error = "stop"))
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


test_that("Processing GO species YAML", {
  expected <- tibble::tribble(
    ~taxon_id, ~full_name, ~common_name_goc, ~code_uniprot, ~taxonomic_group, ~group,
    "NCBITaxon:4932", "Saccharomyces cerevisiae", "budding yeast", "YEAST", "Fungi", "SGD",
    "NCBITaxon:9606", "Homo sapiens", "human", "HUMAN", "Mammalia", "UniProt",
    "NCBITaxon:9999", "Missing code", "missing", NA_character_, "Test group", "UniProt"
  )

  spec <- readr::read_lines("../test_data/goex_test.yaml") |>
    parse_go_species_yaml()

  expect_s3_class(spec, "tbl")
  expect_equal(spec, expected)
})


test_that("Expected behaviour from a non-responsive server", {
  species <- "YEAST-mod"
  httr2::with_mocked_responses(
    mock = mocked_500,
    code = {
      test_unresponsive_server(fetch_go_terms, use_cache = FALSE)
      test_unresponsive_server(fetch_go_species)
      test_unresponsive_server(fetch_go_from_go, species = species, use_cache = FALSE)
      test_unresponsive_server(fetch_go_from_bm, dataset = "scerevisiae_gene_ensembl", use_cache = FALSE)
      test_unresponsive_server(fetch_go_genes_bm, dataset = "scerevisiae_gene_ensembl", use_cache = FALSE)
    })
})


test_that("Expected all_sources behaviour from fetch_go_species", {
  goex <- readr::read_file("../test_data/goex_test.yaml")

  httr2::with_mocked_responses(
    mock = function(req) {
      httr2::response(status_code = 200, body = charToRaw(goex))
    },
    code = {
      spec_all <- fetch_go_species(on_error = "stop", all_sources = TRUE)
      spec_preferred <- fetch_go_species(on_error = "stop", all_sources = FALSE)
    })

  expect_setequal(
    spec_all$designation,
    c("YEAST-mod", "YEAST-uniprot", "HUMAN-uniprot")
  )
  expect_setequal(
    spec_preferred$designation,
    c("YEAST-mod", "HUMAN-uniprot")
  )
  expect_setequal(
    spec_all$annotation_source[spec_all$species == "Saccharomyces cerevisiae"],
    c("MOD", "UniProt")
  )
  expect_false(any(spec_all$species == "Missing code"))
})


test_that("Expected return from fetch_go_species", {
  expected_selection <- c("DROME-mod", "MOUSE-mod", "RAT-mod", "YEAST-mod", "CAEEL-mod", "HUMAN-uniprot")
  spec <- fetch_go_species(on_error = "ignore")
  if(!is.null(spec)) {
    expect_is(spec, "tbl")
    expect_true(all(expected_selection %in% spec$designation))
  }
})


test_that("GO legacy species mapping is exported", {
  legacy <- get_go_legacy_mapping()

  expect_s3_class(legacy, "tbl")
  expect_equal(names(legacy), c("legacy_name", "designation"))
  expect_true(all(c("sgd", "cgd", "xenbase") %in% legacy$legacy_name))
  expect_true("YEAST-mod" %in% legacy$designation)
})


test_that("Correctly match legacy species name", {
  expect_equal(
    match_go_legacy_species("sgd", on_error = "stop"),
    "YEAST-mod"
  )
})


test_that("Ambiguous legacy species names ask users to choose a designation", {
  expect_error(
    match_go_legacy_species("xenbase", on_error = "stop"),
    "corresponds to multiple current designations"
  )

  expect_warning(
    res <- match_go_legacy_species("xenbase", on_error = "warn"),
    "corresponds to multiple current designations"
  )
  expect_null(res)

  expect_null(
    match_go_legacy_species("xenbase", on_error = "ignore")
  )
})


test_that("fetch_go accepts unambiguous legacy GO species names", {
  observed_species <- NULL

  expected_terms <- tibble::tibble(
    term_id = "GO:0000001",
    term_name = "mitochondrion inheritance"
  )
  expected_mapping <- tibble::tibble(
    gene_symbol = "ABC1",
    gene_id = "GENE001",
    db_id = "S000000001",
    term_id = "GO:0000001",
    evidence = "IDA"
  )

  testthat::local_mocked_bindings(
    fetch_go_species = function(on_error = c("stop", "warn", "ignore"),
                                all_sources = FALSE) {
      tibble::tibble(designation = "YEAST-mod")
    },
    fetch_go_genes_go = function(species, use_cache, on_error) {
      observed_species <<- species
      expected_mapping
    },
    fetch_go_terms = function(use_cache, on_error) {
      expected_terms
    }
  )

  res <- fetch_go(species = "sgd", use_cache = FALSE, on_error = "stop")

  expect_equal(observed_species, "YEAST-mod")
  expect_equal(res$terms, expected_terms)
  expect_equal(res$mapping, expected_mapping)
})


test_that("fetch_go rejects ambiguous legacy GO species names", {
  testthat::local_mocked_bindings(
    fetch_go_species = function(on_error = c("stop", "warn", "ignore"),
                                all_sources = FALSE) {
      tibble::tibble(designation = c("XENLA-mod", "XENTR-mod"))
    }
  )

  expect_error(
    fetch_go(species = "xenbase", use_cache = FALSE, on_error = "stop"),
    "corresponds to multiple current designations"
  )
})


test_that("GO GAF files are parsed into gene mapping", {
  expected <- tibble::tribble(
    ~gene_symbol, ~gene_id, ~db_id, ~term_id, ~evidence,
    "ABC1", "GENE001", "S000000001", "GO:0000001", "IDA",
    "DEF2", "GENE002", "S000000002", "GO:0000002", "IEA"
  )

  old_options <- options(
    GO_ANNOTATION_URL = normalizePath("../test_data/annotations", mustWork = TRUE)
  )
  on.exit(options(old_options), add = TRUE)
  testthat::local_mocked_bindings(
    assert_url_path = function(url_path, on_error = c("stop", "warn", "ignore"), timeout = 30) TRUE
  )

  mapping <- fetch_go_genes_go(
    species = "TEST-mod",
    use_cache = FALSE,
    on_error = "stop"
  )

  expect_s3_class(mapping, "tbl")
  expect_equal(mapping, expected)
})


test_that("GO yeast from GO is correct", {
  species <- "YEAST-mod"

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

  re <- fetch_go(species = species, on_error = "ignore")
  if(!is.null(re)) {
    test_fetched_structure(re)
    test_terms(re$terms, expected_terms)
    test_mapping(re$mapping, expected_mapping, "gene_symbol")
  }
})


test_that("GO yeast from Ensembl is correct", {
  dataset <- "scerevisiae_gene_ensembl"

  expected_terms <- tibble::tribble(
    ~term_id, ~term_name,
    "GO:0000166", "nucleotide binding",
    "GO:0006096", "glycolytic process",
    "GO:0005199", "structural constituent of cell wall",
    "GO:0004365", "glyceraldehyde-3-phosphate dehydrogenase (NAD+) (phosphorylating) activity"
  )

  expected_mapping <- tibble::tribble(
    ~term_id, ~gene_id,
    "GO:0000166", "YNL102W",
    "GO:0006096", "YKL060C",
    "GO:0005199", "YKL096W-A",
    "GO:0004365", "YGR192C"
  )

  re <- fetch_go(dataset = dataset, on_error = "ignore")

  if(!is.null(re)) {
    test_fetched_structure(re)
    test_terms(re$terms, expected_terms)
    test_mapping(re$mapping, expected_mapping, "gene_id")
  }
})
