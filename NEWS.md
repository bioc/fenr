## Version 0.1.5

 - added vignette
 - removeLazyData
 - DESCRIPTION updaates
 - replacing url_exists with RCurl::url.exists
 - added NEWS.md
 
## Version 0.1.6 (2022-08-23)

 - bug fix in fetch_kegg
 - removed `min_count` and `fdr_limit` arguments from `functional_enrichment`; filtering can be done afterwards
 - added a small Shiny app as an example of `fenr`
 - updates to vignette

## Version 0.1.7 (2022-09-01)

 - small fixes
 - link to a separate GitHub Shiny app added
 - added support for WikiPathways
 - improved robustness
 - more tests
 
## Version 0.1.8 (2022-09-06)
 
 - `fetch_reactome` provides two ways of retrieving data, via one downloadable file or via APIs
 
## Version 0.1.9 (2022-09-13)
 
  - Replacing ontologyIndex::get_ontology with a simpler parser
  - Replacing KEGGREST with own simple API parsers
  - Applying BiocStyle to the vignette
  - improving test_functional_enrichment
  
## Version 0.1.10 (2022-10-04)

 - Style changes for BiocCheck
 - Adding more tests
 - Fixing a bug in `parse_kegg_genes`
 
## Version 0.1.11 (2022-10-11)
 
 - Significant speed-up of enrichment by using Rfast::Hash in place of R lists
 - KEGG improvements, recognizing flat file genes with no gene synonym
 - Additional tests for Reactome
 - Minor improvements and fixes

## Version 0.1.13
 
- Added functions `get_term_features` and `get_feature_terms` to access data safely
- HACK: BioPlanet server's SSL certificate expired, so need insecure download.

## Version 0.1.14 (2023-02-02)

- Ditched large and clunky `Rfast` and using native R environments as fast hashes (see https://riptutorial.com/r/example/18339/environments-as-hash-maps)
- A few tweaks and improvements

## Version 0.1.15 (2023-03-02)

- Fixed a bug where there are some features at a term that are not present in the universe (all features). This could happen when the universe was particularly small. Potentially a serious bug.

## Version 0.1.16 (2023-03-16)

- Added a fix to work correctly with integer feature IDs.


## Version 0.1.17 (2023-04-18)

 - New units tests added

## Version 0.99.0 (2023-04-20)

 - Pre-release Bioconductor version.

## Version 0.99.1 (2023-04-25)

 - BioPlanet database vanished from internet and there is no sign of it coming back. Removing all BioPlanet-related code and replacing BioPlanet with GO in the vignette and examples (this, alas, makes it longer to check).
 - OK, it is back, but I keep GO examples and vignettes.
 - Minor improvements to documentation.
 
## Version 0.99.2 (2023-05-24)
 
 - BioPlanet's tripod.nih.gov SSL certificate seems to be fixed, so reversing to the original read_csv code.
 
## Version 0.99.3
 
 - Continuing issues with access to BioPlanet. `fetch_bp` example is now marked `donotrun` and testing `fetch_bp` is removed to ensure smooth build and check even when BioPlanet server is down.
 
## Version 0.99.4
 
 - Taking `BiocCheck` new warnings into account: adding `@return` to data roxygens.
 
## Version 0.99.5

 - Major overhaul following comments from Bioconductor's reviewer.

## Version 0.99.6

### In response to reviewer's comments

 - The wording in the vignette was adjusted to more clearly convey the purpose of the package to users
 - Rewritten the description in DESCRIPTION file to clearly convey the purpose of the package to users
 
### BioPlanet seems defunct 
 
 - Removed BioPlanet for good, as their webpage is continuously down and the maintainer is not responding

### Minor adjustments to speed up building and testing

 - Removed KEGG from interactive example to speed up vignette building (GO and Reactome are sufficient for a simple example)
 - Replaced yeast with simpler organisms in Wiki and KEGG tests to speed up testing
 - Replaced yeast with simpler organisms in Wiki and KEGG examples to speed up checking
 
## Version 0.99.7

 - Minor changes to prepare for Bioconductor release
 - Reverting temporarily to *readr* version 1 to circumvent a *vroom* 1.6.4 bug

## Version 1.1.1

 - First update after Bioconductor release
 - Implemented changes to prevent the package from build/check fail, if one of the remote servers is not responding
 - Moved from `httr` to `httr2` 
 - Tests and examples now generate warnings in case of server failure
 - Added tests for behaviour in case of a non-responsive server
 - Extended test coverage to 100%, except for the interactive example
 
## Version 1.1.2

 - Bug fixes, examples need on_error = "warn"

