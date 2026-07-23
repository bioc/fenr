#' Example set of background genes.
#'
#' A set of gene names for proteins that could be detected in a typical
#' proteomics experiment on yeast samples.
#'
#' @usage data(exmpl_all)
#' @format A character vector with 6985 elements.
#' @return A vector with background gene names.
"exmpl_all"


#' Example set of selected genes.
#'
#' A set of gene names manually selected to illustrate functional enrichment.
#'
#' @usage data(exmpl_sel)
#' @format A character vector with 21 elements.
#' @return A vector with selected gene names.
"exmpl_sel"


#' Differential expression results for yeast RNA-seq.
#'
#' A subset of 6 + 6 replicates was selected from data set reported in
#' https://doi.org/10.1093/bioinformatics/btv425
#'
#' @usage data(yeast_de)
#' @format A tibble with 5 columns
#' @return Results for differential expression for yeast RNA-seq.
"yeast_de"


#' GO-terms data downloaded for the vignette.
#'
#' Downloaded using \code{go <- fetch_go(species = "YEAST-mod")}
#'
#' @usage data(go)
#' @format A list of two tibbles
#' @return Contains GO-term descriptions and gene mapping.
"go"


#' GO species
#'
#' Downloaded using \code{go_species <- fetch_go_species()}
#'
#' @usage data(go_species)
#' @format A tibble
#' @return Contains species available through Gene Ontology
"go_species"
