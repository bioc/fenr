#' Create term data for interactive example
#'
#' @param de Differential expression results, use \code{yeast_de} data attached
#'   to this package.
#'
#' @return A list of objects containing functional terms for GO and Reactome.
#' @export
#' @examples
#' data(yeast_de)
#' term_data <- fetch_terms_for_example(yeast_de)
fetch_terms_for_example <- function(de) {
  # All gene background
  all_genes <- de$gene_id

  # load GO terms
  message("Fetching GO data\n")
  go <- fetch_go(species = "sgd")
  go_data <- prepare_for_enrichment(go$terms, go$mapping, all_genes, feature_name = "gene_synonym")

  # load Reactome pathways
  message("Fetching Reactome data\n")
  re <- fetch_reactome("Saccharomyces cerevisiae", on_error = "warn")
  re_data <- prepare_for_enrichment(re$terms, re$mapping, all_genes, feature_name = "gene_id")

  # Put all functional term data in one structure; Shiny app will access
  # individual ontologies from this list
  list(
    go = go_data,
    re = re_data
  )
}


#' Volcano plot
#'
#' @param d Tibble with x, y and FDR
#' @param fdr_limit FDR limit below which point will be plotted black
#'
#' @return A ggplot object
#' @noRd
plot_volcano <- function(d, fdr_limit = 0.05) {
  # Binding variables from non-standard evaluation locally
  x <- y <- FDR <- NULL

  sres <- d |>
    dplyr::filter(FDR <= fdr_limit)
  g <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 18)
    ) +
    ggplot2::geom_point(size = 0.2) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70") +
    ggplot2::labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.03)))
  if (nrow(sres) > 0) {
    g <- g + ggplot2::geom_hline(yintercept = -log10(max(sres$PValue)), colour = "red",
                        linetype = "dashed", alpha = 0.2)
  }
  g
}

#' MA plot
#'
#' @param d Tibble with x, y and FDR
#' @param fdr_limit FDR limit below which point will be plotted black
#'
#' @return A ggplot object
#' @noRd
plot_ma <- function(d, fdr_limit = 0.05) {
  # Binding variables from non-standard evaluation locally
  x <- y <- NULL

  ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 18)
    ) +
    ggplot2::geom_point(data = d[d$FDR > fdr_limit, ], size = 0.1, colour = "grey50") +
    ggplot2::geom_point(data = d[d$FDR <= fdr_limit, ], size = 0.2, colour = "black") +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::labs(x = expression(log[2]~CPM), y = expression(log[2]~FC))
}



#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param input App input
#'
#' @return A tibble with x, y
#' @noRd
get_xy_data <- function(de, input) {
  logFC <- PValue <- logCPM <- NULL
  if (input$plot_type == "volcano") {
    xy_data <- de |>
      dplyr::mutate(x = logFC, y = -log10(PValue))
  } else if (input$plot_type == "ma") {
    xy_data <- de  |>
      dplyr::mutate(x = logCPM, y = logFC)
  }
  xy_data
}

#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param term_data All functional term data
#' @param input App input
#' @param max_points Maximum number of points allowed to select
#'
#' @return A tibble with functional enrichment results
#' @noRd
enrichment_table <- function(de, term_data, input, max_points = 10000) {
  # Binding variables from non-standard evaluation locally
  p_adjust <- n_with_sel <- NULL

  xy_data <- get_xy_data(de, input)
  terms <- term_data[[input$ontology]]
  gene2symbol <- rlang::set_names(de$gene_symbol, de$gene_id)
  sel <- NULL
  fe <- NULL
  if (!is.null(input$plot_brush)) {
    brushed <- shiny::brushedPoints(xy_data, input$plot_brush)
    sel <- brushed$gene_id
    n <- length(sel)
    if (n > 0 && n <= max_points) {
      fe <- functional_enrichment(de$gene_id, sel, terms, feat2name = gene2symbol)
      if(!is.null(fe))
        fe <- dplyr::filter(fe, p_adjust < 0.05 & n_with_sel > 2)
    } else if (n > 0) {
      fe <- data.frame(Error = paste("only", max_points, "points can be selected."))
    }
  }
  fe
}

#' Helper function for Shiny app
#'
#' @param de Differential expression results
#' @param input App input
#'
#' @return Data for volcano or ma plot
#' @noRd
main_plot <- function(de, input) {
  xy_data <- get_xy_data(de, input)
  if (input$plot_type == "volcano") {
    plot_volcano(xy_data)
  } else if (input$plot_type == "ma") {
    plot_ma(xy_data)
  }
}


# Note for Bioconductor: Using "donttest" is problematic here. When running
# `devtools::check()`, it stalls at "checking examples with --run-donttest ...".
# I suspect it's attempting to launch the Shiny app, which fails. As a result,
# contrary to `BiocCheck()` guidelines, we've opted for "dontrun" in this instance.


#' Small Shiny app serving as example for fast enrichment
#'
#' @param de Differential expression results, \code{yeast_de} attached to this
#'   package can be used.
#' @param term_data A list of \code{fenr_terms} objects containing functional
#'   data for various ontologies. \code{fetch_terms_for_example} can be used to
#'   create this object.
#'
#' @return An interactive Shiny app
#' @importFrom assertthat assert_that
#' @export
#' @examples
#' data(yeast_de)
#' term_data <- fetch_terms_for_example(yeast_de)
#' if(interactive()) {
#'   enrichment_interactive(yeast_de, term_data)
#' }
enrichment_interactive <- function(de, term_data) {
  assert_that(is.data.frame(de))
  assert_that(is(term_data, "list"))

  req_cols <- c("gene_id", "gene_symbol", "logFC", "logCPM", "PValue", "FDR")
  if(!all(req_cols %in% colnames(de)))
    stop(paste0("The data frame needs to contain the following columns\n", paste(req_cols, collapse = ", ")))


  ui <- function() {
    shiny::shinyUI(shiny::fluidPage(
      shiny::tags$style("table{font-size: 11px; background-color: #EAF5FF}"),
      shiny::titlePanel("fenr example"),
      shiny::p("Select a group of genes in the plot to see their functional enrichment."),

      shiny::fluidRow(
        shiny::column(5,
          shiny::radioButtons(
            inputId = "plot_type",
            label = "Plot type:",
            choices = c("Volcano" = "volcano", "MA" = "ma"),
            inline = TRUE
          ),
          shiny::plotOutput(
            outputId = "main_plot",
            height = "480px",
            width = "100%",
            brush = "plot_brush",
            hover = "plot_hover"
          )
        ),
        shiny::column(7,
          shiny::radioButtons(
            inputId = "ontology",
            label = "Ontology:",
            choices = c("GO" = "go", "Reactome" = "re"),
            inline = TRUE
          ),
          shiny::div(style = "height: 480px; overflow-y: scroll", shiny::tableOutput("enrichment")),
        )
      )
    ))
  }

  server <- function(input, output, session) {
    # Prevents RStudio from crashing when Shiny window closed manually
    session$onSessionEnded(function() {
      shiny::stopApp()
    })

    output$enrichment <- shiny::renderTable({
      enrichment_table(de, term_data, input)
    })

    output$main_plot <- shiny::renderPlot({
      main_plot(de, input)
    })
  }

  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}
