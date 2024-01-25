fpath <- "https://raw.githubusercontent.com/bartongroup/fenr/main/DESCRIPTION"
rname <- "test"
expected_content <- readLines(fpath)

cleanup_test <- function(bfc) {
  hits <- BiocFileCache::bfcquery(bfc, rname, field = "rname")
  if(nrow(hits) > 0)
    BiocFileCache::bfcremove(bfc, hits$rid)
}


test_that("Cache location exists", {
  cache <- cache_location()
  expect_true(dir.exists(cache))
})


test_that("No caching test", {
  re <- cached_url_path(rname, fpath, use_cache = FALSE)
  expect_equal(re, fpath)
})


test_that("Adding new entry to the cache", {
  cache <- cache_location()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  cleanup_test(bfc)

  re <- cached_url_path(rname, fpath, use_cache = TRUE)
  expect_true(file.exists(re))
  re_content <- readLines(re)
  expect_equal(re_content, expected_content)
})


test_that("Retrieving an entry to the cache", {
  cache <- cache_location()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  cleanup_test(bfc)
  lpath <- BiocFileCache::bfcadd(bfc, rname, fpath)

  re <- cached_url_path(rname, fpath, use_cache = TRUE)
  expect_true(file.exists(re))
  expect_equal(re, lpath)

  re_content <- readLines(re)
  expect_equal(re_content, expected_content)
})


test_that("Multiple entries behaviour", {
  cache <- cache_location()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  cleanup_test(bfc)
  lpath1 <- BiocFileCache::bfcadd(bfc, rname, fpath)
  lpath2 <- BiocFileCache::bfcadd(bfc, rname, fpath)

  re <- cached_url_path(rname, fpath, use_cache = TRUE)
  expect_true(file.exists(re))

  re_content <- readLines(re)
  expect_equal(re_content, expected_content)
})


test_that("Clearing cache", {
  remove_cache(ask = FALSE)

  cache <- cache_location()
  bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
  cnt <- BiocFileCache::bfccount(bfc)

  expect_equal(cnt, 0)
})
