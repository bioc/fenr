# Fast functional enrichment

Maintainer: Marek Gierlinski (<M.Gierlinski@dundee.ac.uk>)

This R package provides a fast and efficient method for functional enrichment analysis, optimized for speed and designed for use in interactive applications, such as *Shiny* apps.


## Quick example

To run a *Shiny* app demonstration of `fenr` directly from GitHub, enter the following command in your R console:

```
shiny::runGitHub("bartongroup/fenr-shiny-example")
```

## Installation

`fenr` can be installed using

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("fenr")
```

## Usage

The initial step involves downloading functional term data. `fenr` supports data downloads from *Gene Ontology*, *Reactome*, *KEGG*, and *WikiPathways*. Custom ontologies can also be used, provided they are converted into an appropriate format (refer to the `prepare_for_enrichment` function for more information). The command below downloads functional terms and gene mapping from Gene Ontology (GO) for yeast:

```
go <- fetch_go(species = "YEAST-mod")
```

This command returns a list with two tibbles containing term information (`term_id` and `term_name`) and gene-term mapping (`term_id` and `gene_symbol`). We convert this data into an object suitable for fast functional enrichment. `exmpl_all` is an example of gene background provided by the package, which contains a vector with gene symbols related to all detections in an experiment.

```
data(exmpl_all, exmpl_sel)
go_terms <- prepare_for_enrichment(go$terms, go$mapping, exmpl_all, feature_name = "gene_symbol")
```

The `go_terms` object is a data structure containing all mappings in a quickly accessible form. From this point on, you can use go_terms to perform multiple functional enrichments on various gene selections. For example, if `exmpl_all` is a vector with all background gene symbols and `exmpl_sel` is a vector with genes of interest (both provided by the package), you can perform functional enrichment analysis using:

```
enr <- functional_enrichment(exmpl_all, exmpl_sel, go_terms)
```

The result is a tibble:

```
# A tibble: 51 × 10
  term_id    term_name                                N_with n_with_sel n_expect enrichment odds_ratio ids    p_value p_adjust
   <chr>      <chr>                                     <int>      <int>    <dbl>      <dbl>      <dbl> <chr>    <dbl>    <dbl>
 1 GO:1905356 regulation of snRNA pseudouridine synth…      2          2     0.01      333        Inf   TOR2… 8.61e- 6 5.49e- 5
 2 GO:0031929 TOR signaling                                19         18     0.06      315      41800   TOR2… 0        0       
 3 GO:0031931 TORC1 complex                                11         10     0.03      302       6330   TOR2… 0        0       
 4 GO:0001558 regulation of cell growth                     9          5     0.03      185        544   KOG1… 1.84e-11 2.34e-10
 5 GO:0031932 TORC2 complex                                15          7     0.05      155        435   TOR2… 4.66e-15 7.93e-14
 6 GO:0016242 negative regulation of macroautophagy         9          3     0.03      111        193   KSP1… 1.95e- 6 1.42e- 5
 7 GO:0043666 regulation of phosphoprotein phosphatas…     13          4     0.04      102        182   SAP1… 4.24e- 8 4.33e- 7
 8 GO:0030950 establishment or maintenance of actin c…     15          4     0.05       88.7      149   TOR2… 8.07e- 8 6.86e- 7
 9 GO:0031930 mitochondria-nucleus signaling pathway        9          2     0.03       73.9      105   TOR1… 3.06e- 4 1.3 e- 3
10 GO:0010507 negative regulation of autophagy             12          2     0.04       55.4       73.2 TOR2… 5.58e- 4 2.03e- 3
# ℹ 41 more rows
```

The columns are as follows

 - `N_with`: The number of features (genes) associated with this term in the background of all genes.
 - `n_with_sel`: The number of features associated with this term in the selection.
 - `n_expect`: The expected number of features associated with this term under the null hypothesis (terms are randomly distributed).
 - `enrichment`: The ratio of observed to expected.
 - `odds_ratio`: The effect size, represented by the odds ratio from the contingency table.
 - `ids`: The identifiers of features with the term in the selection.
 - `p_value`: The raw p-value from the hypergeometric distribution.
 - `p_adjust`: The p-value adjusted for multiple tests using the Benjamini-Hochberg approach.

# Interactive Example

A small Shiny app is included in the package to demonstrate the usage of `fenr` in an interactive environment. All time-consuming data loading and preparation tasks are performed before the app is launched.

```
data(yeast_de)
term_data <- fetch_terms_for_example(yeast_de)
```
 
`yeast_de` is the result of differential expression (using `edgeR`) on a subset of 6+6 replicates from [Gierlinski et al. (2015)](https://academic.oup.com/bioinformatics/article/31/22/3625/240923).

The function `fetch_terms_for_example` uses `fetch_*` functions from `fenr` to download and process data from *GO*, *Reactome* and *KEGG*.  You can view the step-by-step process by examining the function code on [GitHub](https://github.com/bartongroup/fenr/blob/main/R/iteractive_example.R). The object `term_data` is a named list of `fenr_terms` objects, one for each ontology.

After completing the slow tasks, you can start the Shiny app by running:

```
enrichment_interactive(yeast_de, term_data)
```

To quickly see how `fenr` works an example can be loaded directly from GitHub:

```
shiny::runGitHub("bartongroup/fenr-shiny-example")
```
