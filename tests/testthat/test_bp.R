expected_mapping <- tibble::tribble(
  ~term_id, ~gene_symbol,
  "bioplanet_1025", "CDK1",
  "bioplanet_120", "IL1A",
  "bioplanet_1755", "RPL6",
  "bioplanet_1121", "POLR1A"
)

test_that("Bioplanet mapping makes sense", {
  bp <- fetch_bp(on_error = "warn")
  if(!is.null(bp)) {
    expect_is(bp, "list")
    expect_setequal(names(bp), c("terms", "mapping"))
    mapping <- bp$mapping

    merged <- expected_mapping |>
      dplyr::left_join(mapping, by = c("term_id", "gene_symbol")) |>
      tidyr::drop_na()

    expect_equal(nrow(expected_mapping), nrow(merged))
  }
})


test_that("Expected behaviour from a non-responsive server", {
  httr2::with_mocked_responses(
    mock = mocked_500,
    code = {
      test_unresponsive_server(fetch_bp, use_cache = FALSE)
    })
})
