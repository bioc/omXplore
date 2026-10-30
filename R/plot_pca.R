#' @title PCA plots
#' @importFrom stats na.omit
#'
#' @param res.pca The result of the function `FactoMineR::PCA()`
#' @param qdata A data.frame() of quantitative data
#' @param group A vector with the name of samples
#' @param var.scaling A boolean indicating whether to scale the data or not
#' @param ncp See `FactoMineR::PCA()`
#' @param chosen.axes See the parameter 'axes' of the function
#' `factoextra::fviz_pca_var()`
#'
#' @name ds-pca
#'
#' @examples
#' if (interactive()) {
#'     data(vdata)
#'     obj <- vdata[[1]]
#'     res.pca <- wrapper_pca(qdata = SummarizedExperiment::assay(obj), group = get_group(obj))
#'     plotPCA_Eigen(res.pca)
#'     plotPCA_Var(res.pca)
#'     plotPCA_Eigen_hc(res.pca)
#'     plotPCA_Ind(res.pca)
#' }
#'
NULL




#' @export
#' @return The result of the function FactoMineR::PCA()
#' @rdname ds-pca
#' @import FactoMineR
#'
wrapper_pca <- function(
        qdata,
        group,
        var.scaling = TRUE,
        ncp = NULL,
        approach = "FM",
        gramschmidt = TRUE) {
    if (missing(qdata)) {
        stop("'qdata' is missing.")
    }

    stopifnot(inherits(qdata, "matrix"))


    if (is.null(var.scaling)) {
        var.scaling <- TRUE
    }

    res.pca <- NULL

    # if (length(which(is.na(obj@qdata))) > 0) {
    if (is.null(ncp)) {
        nmax <- 12
        y <- qdata
        nprot <- dim(y)[1]
        n <- dim(y)[2] # If too big, take the number of conditions.

        if (n > nmax) {
            n <- length(unique(group))
        }

        ncp <- min(n, nmax)
    }


    res.pca <- my_PCA(
        X = qdata,
        scale.unit = var.scaling,
        ncp = ncp,
        graph = FALSE,
        approach = approach,
        gramschmidt = gramschmidt
    )

    return(res.pca)
}




#' @export
#' @import plotly
#'
#' @rdname ds-pca
#' @return A plot
#'
plotPCA_Eigen <- function(res.pca) {
    stopifnot(!is.null(res.pca))

    pc_names <- rownames(res.pca$eig)
    
    df <- data.frame(
        PC = factor(pc_names, levels = pc_names),
        var = res.pca$eig[, 2],
        cumvar = res.pca$eig[, 3],
        eig = res.pca$eig[, 1]
    )
    
    
    p <- plotly::plot_ly(df, x = ~PC) |>
        plotly::add_bars(
            y = ~var,
            name = "% of variances",
            marker = list(color = "rgba(100,150,200,0.8)")
        ) |>
        plotly::add_lines(
            y = ~cumvar,
            mode = "lines+markers",
            name = "Cumulative % of variances",
            line = list(color = "darkblue", width = 2),
            marker = list(color = "darkblue", size = 6)
        ) |>
        plotly::layout(
            xaxis = list(
                title = "Principal Components",
                categoryorder = "array",
                categoryarray = pc_names),
            yaxis = list(
                title = "% of variances",
                range = c(0, 100),
                ticksuffix = "%"),
            margin = list(t = 30),
            legend = list(
                orientation = "h",
                x = 0,
                y = -0.15
            )
        )
    
    return(p)
}





#' @return A plot
#'
#' @rdname ds-pca
#' @export
#' @import factoextra
#'
plotPCA_Var <- function(res.pca, chosen.axes = c(1, 2)) {
    # plot.PCA(res.pca, choix="var", axes = chosen.axes,
    # title="Sample factor map (PCA)")
    # Colorer en fonction du cos2: qualite de representation
    if (is.null(res.pca)) {
        return(NULL)
    }
    factoextra::fviz_pca_var(
        res.pca,
        axes = chosen.axes,
        col.var = "cos2",
        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
        repel = TRUE # Evite le chevauchement de texte
    )
}




#' @export
#' @rdname ds-pca
#' @return A plot
#' @import factoextra
#'
plotPCA_Ind <- function(res.pca, chosen.axes = c(1, 2)) {
    if (is.null(res.pca)) {
        return(NULL)
    }

    factoextra::fviz_pca_ind(res.pca,
        axes = chosen.axes,
        geom = "point"
    )
}





#' @import plotly
#' @rdname ds-pca
#' @return A plot
#' @export
#'
plotPCA_Eigen_hc <- function(res.pca) {
    if (is.null(res.pca)) {
        return(NULL)
    }
    
    pc_names <- rownames(res.pca$eig)
    
    df <- data.frame(
        PC = factor(pc_names, levels = pc_names),
        var = res.pca$eig[, 2],
        cumvar = res.pca$eig[, 3],
        eig = res.pca$eig[, 1]
    )
    
    
    p <- plotly::plot_ly(df, x = ~PC) |>
        plotly::add_bars(
            y = ~var,
            name = "% of variances",
            marker = list(color = "rgba(100,150,200,0.8)")
        ) |>
        plotly::add_lines(
            y = ~cumvar,
            mode = "lines+markers",
            name = "Cumulative % of variances",
            line = list(color = "darkblue", width = 2),
            marker = list(color = "darkblue", size = 6)
        ) |>
        plotly::layout(
            xaxis = list(
                title = "Principal Components",
                categoryorder = "array",
                categoryarray = pc_names),
            yaxis = list(
                title = "% of variances",
                range = c(0, 100),
                ticksuffix = "%"),
            margin = list(t = 30),
            legend = list(
                orientation = "h",
                x = 0,
                y = -0.15
            )
        )
    
    return(p)
}
