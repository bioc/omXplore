#' @title Bar plot of missing values per lines using plotly.
#'
#' @description
#'
#' This method plots a bar plot which represents the distribution of the
#' number of missing values (NA) per lines (ie proteins).
#'
#' * `wrapper_pca()`
#' * `plotPCA_Eigen_hc()`: plots the eigen values of PCA with the plotly
#'    library
#' * `plotPCA_Eigen()`: plots the eigen values of PCA
#' * `plotPCA_Var()`
#' * `plotPCA_Ind()`
#'
#' @param id A `character(1)` which is the id of the shiny module.
#' @param dataIn An instance of the class `MultiAssayExperiment`.
#' @param i An integer which is the index of the assay in the param obj
#' @param var.scaling The dimensions to plot
#' @param ncp A `integer(1)` which represents the umber of dimensions kept in
#' the results.
#' @param res.pca Result of FactoMineR::PCA
#' @param chosen.axes The dimensions to plot
#'
#' @author Samuel Wieczorek, Enora Fremy
#'
#' @name ds-pca
#'
#'
#' @author Samuel Wieczorek, Enora Fremy
#'
#' @examples
#' if (interactive()) {
#'     data(vdata)
#'     library(shiny)
#'     library(dplyr)
#'     library(plotly)
#'     # Replace missing values for the example
#'     sel <- is.na(SummarizedExperiment::assay(vdata, 1))
#'     SummarizedExperiment::assay(vdata[[1]])[sel] <- 0
#'     SummarizedExperiment::assay(vdata[[1]])[1, 1] <- NA
#'     omXplore_pca(vdata, 1)
#' }
#'
NULL



#' @import shiny
#' @importFrom shinyjs useShinyjs hidden toggle
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly renderPlotly
#' @importFrom shinyjs useShinyjs hidden toggle
#'
#' @rdname ds-pca
#' @export
#' @return NA
#'
omXplore_pca_ui <- function(id) {
    ns <- NS(id)
    tagList(
        shinyjs::useShinyjs(),
        shinyjs::hidden(div(
            id = ns("badFormatMsg"),
            h3(globals()$bad_format_txt)
        )),
        uiOutput(ns("WarningNA_PCA")),
        uiOutput(ns("pcaOptions")),
        tags$hr(),
        uiOutput(ns("pcaPlots"))
    )
}


#' @import shiny
#' @importFrom shinyjs useShinyjs hidden toggle
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotly renderPlotly
#' @importFrom shinyjs useShinyjs hidden toggle
#' @importFrom SummarizedExperiment assay
#'
#' @rdname ds-pca
#'
#' @import factoextra
#'
#' @export
#' @return NA
#'
omXplore_pca_server <- function(
        id,
        dataIn,
        i) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        rv.pca <- reactiveValues(
            data = NULL,
            PCA_axes = NULL,
            res.pca = NULL,
            PCA_varScale = TRUE,
            gramschmidt_PCA = TRUE,
            approach_PCA = "FM",
            i = NULL
        )

        observeEvent(dataIn(), {
            is.mae <- inherits(dataIn(), "MultiAssayExperiment")
            stopifnot(is.mae)
            req(i())
            if (i() %in% names(dataIn())){
                rv.pca$i <- i()
            } else {
                rv.pca$i <- names(dataIn())[length(dataIn())]
            }
            
            rv.pca$data <- SummarizedExperiment::assay(dataIn(), rv.pca$i)
            shinyjs::toggle("badFormatMsg", condition = !is.mae)
        })

        output$WarningNA_PCA <- renderUI({
            # req(rv.pca$data)
            # req(length(which(is.na(rv.pca$data))) > 0)

            rule <- !is.null(rv.pca$data) && length(which(is.na(rv.pca$data))) > 0
            req(rule)
            if (rv.pca$approach_PCA == "FM") {
                tagList(
                    tags$p(
                        style = "color:red;font-size: 16px",
                        "Warning: PCA cannot be computed with FactoMineR as your dataset contains missing values. 
                        Please impute the them first or use NIPALS."
                    )
                )
            }
        })
        
        output$PCA_varScale_UI <- renderUI({
            req(rv.pca$approach_PCA == "FM")
            
            checkboxInput(ns("PCA_varScale"), 
                          "VarScale in FactoMineR", 
                          value = TRUE)
        })
            
        output$gramschmidt_PCA_UI <- renderUI({
            req(rv.pca$approach_PCA == "NIPALS")
            
            checkboxInput(ns("gramschmidt_PCA"), 
                          "gramschmidt in Nipals", 
                          value = TRUE)
        })


        output$pcaOptions <- renderUI({
            req(rv.pca$data)
            #print(length(rv.pca$data))
            # req(length(which(is.na(rv.pca$data))) == 0)
            tagList(
                tags$div(style = "display: flex; gap: 8px;",
                     tags$div(style = "white-space: nowrap; margin-right: 25px;",
                         selectInput(ns("approach_PCA"),
                                     label = "Approach used for PCA",
                                     choices = c("FactoMineR" = "FM", "Nipals" = "NIPALS"),
                                     width = "150px"
                         )
                     ),
                    tags$div(
                        numericInput(ns("pca_axe1"), "Dimension 1",
                            min = 1,
                            max = Compute_PCA_dim(),
                            value = 1,
                            width = "100px"
                        )
                    ),
                    tags$div(
                        numericInput(ns("pca_axe2"), "Dimension 2",
                            min = 1,
                            max = Compute_PCA_dim(),
                            value = 2,
                            width = "100px"
                        )
                    ),
                    div(style = "margin-top: 10px; margin-left: 15px;",
                        uiOutput(ns("PCA_varScale_UI")),
                        uiOutput(ns("gramschmidt_PCA_UI"))
                    )
                )
            )
        })

        observeEvent(c(input$pca_axe1, input$pca_axe2), {
            rv.pca$PCA_axes <- c(input$pca_axe1, input$pca_axe2)
        })
        observeEvent(req(input$PCA_varScale), {
            rv.pca$PCA_varScale <- input$PCA_varScale
        })
        observeEvent(input$gramschmidt_PCA, {
            rv.pca$gramschmidt_PCA <- input$gramschmidt_PCA
        })
        observeEvent(req(input$approach_PCA), {
            print("ava")
            print(rv.pca$approach_PCA)
            rv.pca$approach_PCA <- input$approach_PCA
            print("apr")
            print(rv.pca$approach_PCA)
        })
        observeEvent(input$PCA_varScale, {
            rv.pca$PCA_varScale <- input$PCA_varScale
        })

        observe({
            rule1 <- rv.pca$approach_PCA == "FM" && length(which(is.na(rv.pca$data))) == 0
            rule2 <- rv.pca$approach_PCA == "NIPALS"
            req(rule1 || rule2)

            rv.pca$res.pca <- wrapper_pca(
                qdata = SummarizedExperiment::assay(dataIn(), rv.pca$i),
                group = get_group(dataIn()),
                var.scaling = rv.pca$PCA_varScale,
                ncp = Compute_PCA_dim(),
                approach = rv.pca$approach_PCA,
                gramschmidt = rv.pca$gramschmidt_PCA
            )
        })

        output$pcaPlots <- renderUI({
            req(rv.pca$data)
            req(rv.pca$res.pca$var$coord)

            rule1 <- rv.pca$approach_PCA == "FM" && length(which(is.na(rv.pca$data))) == 0
            rule2 <- rv.pca$approach_PCA == "NIPALS"
            req(rule1 || rule2)

            tagList(
                plotOutput(ns("pcaPlotVar")),
                plotOutput(ns("pcaPlotInd")),
                div(style = "background-color: white;",
                    formatDT_ui(ns("PCAvarCoord"))),
                plotly::plotlyOutput(ns("pcaPlotEigen"))
            )
        })

        observe({
            df <- as.data.frame(rv.pca$res.pca$var$coord)
            formatDT_server("PCAvarCoord",
                data = reactive({
                    round(df, digits = 2)
                }),
                showRownames = TRUE
            )
        })

        output$pcaPlotVar <- renderPlot({
            req(c(rv.pca$PCA_axes, rv.pca$res.pca))

            withProgress(message = "Making plot", value = 100, {
                factoextra::fviz_pca_var(rv.pca$res.pca,
                    axes = rv.pca$PCA_axes,
                    col.var = "cos2",
                    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                    repel = TRUE
                )
            })
        })

        output$pcaPlotInd <- renderPlot({
            req(c(rv.pca$PCA_axes, rv.pca$res.pca))

            withProgress(message = "Making plot", value = 100, {
                factoextra::fviz_pca_ind(rv.pca$res.pca,
                    axes = rv.pca$PCA_axes,
                    geom = "point"
                )
            })
        })


        output$pcaPlotEigen <- plotly::renderPlotly({
            req(rv.pca$res.pca)

            withProgress(message = "Making plot", value = 100, {
                plotPCA_Eigen(rv.pca$res.pca)
            })
        })


        Compute_PCA_dim <- reactive({
            req(rv.pca$data)
            nmax <- 12 # ncp should not be greater than...
            # for info, ncp = number of components or dimensions in PCA results

            y <- rv.pca$data
            nprot <- dim(y)[1]
            n <- dim(y)[2] # If too big, take the number of conditions.

            if (n > nmax) {
                n <- length(unique(get_group(dataIn())))
            }

            ncp <- min(n, nmax)
            ncp
        })
    })
}


#' @export
#' @rdname ds-pca
#' @return A shiny app
#'
omXplore_pca <- function(dataIn, i) {
    stopifnot(inherits(dataIn, "MultiAssayExperiment"))

    ui <- omXplore_pca_ui("plot")

    server <- function(input, output, session) {
        omXplore_pca_server("plot",
            dataIn = reactive({
                dataIn
            }),
            i = reactive({
                i
            })
        )
    }
    shinyApp(ui = ui, server = server)
}
