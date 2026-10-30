#' @title Display a CC
#'
#' @param cc A Connected Component (a list)
#' @param g An instance of a graph
#' @param layout A `character(1)` which i the layout used in visNetwork.
#' Default value is 'layout_with_fr'
#' @param df A data.frame()
#'
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' if (interactive()) {
#'     data(sub_R25)
#'     se <- sub_R25[[1]]
#'     g <- buildGraph(get_cc(se)[[1]])
#'     display.CC.visNet(g)
#' }
#'
#' @name pep_prot_CC
#'
#' @import plotly
#' @import visNetwork
#'
NULL


#' @rdname pep_prot_CC
#' @return A list
#' @export
buildGraph <- function(cc) {
    # cc is the sparse sub-matrix of one connected component
    subX <- as.matrix(cc)
    
    nb.pep  <- nrow(subX)
    nb.prot <- ncol(subX)
    nb.total <- nb.prot + nb.pep
    
    pep.deg <- rowSums(subX != 0)   # number of proteins per peptide
    
    edge.list <- as.data.frame(which(subX != 0, arr.ind = TRUE))
    
    # Peptide groups: shared if linked to >1 protein, specific if exactly 1
    def.grp <- c(
        ifelse(pep.deg == 1, "spec.peptide", "shared.peptide"),
        rep("protein", nb.prot)
    )
    
    nodes <- data.frame(
        id = seq_len(nb.total),
        group = def.grp,
        label = c(rownames(subX), colnames(subX)),
        title = paste0("<p>", seq_len(nb.total), "<br>Tooltip !</p>"),
        size = c(rep(10, nb.pep), rep(20, nb.prot)),
        stringsAsFactors = FALSE
    )
    
    edges <- data.frame(
        from = edge.list$row,
        to = edge.list$col + nb.pep,
        stringsAsFactors = FALSE
    )
    
    list(nodes = nodes, edges = edges)
}





#' @rdname pep_prot_CC
#' @export
#' @return A plot
#'
display.CC.visNet <- function(
        g = NULL,
        layout = "layout_with_fr") {
    if (is.null(g)) {
        message("g is NULL. Abort...")
        return(NULL)
    }

    col.prot <- "#ECB57C"
    col.spec <- "#5CA3F7"
    col.shared <- "#0EA513"

    visNetwork::visNetwork(g$nodes, g$edges, width = "100%", height = "100%") |>
        visNetwork::visNodes(shape = "dot") |> # square for all nodes
        visNetwork::visGroups(
            groupname = "spec.peptide",
            color = col.spec
        ) |> # darkblue for group "A"
        visNetwork::visGroups(
            groupname = "shared.peptide",
            color = col.shared
        ) |> # darkblue for group "A"
        visNetwork::visGroups(
            groupname = "protein",
            color = col.prot, shape = "dot"
        ) |>
        visNetwork::visOptions(highlightNearest = FALSE) |>
        # visLegend()
        # visPhysics(stabilization = FALSE)|>
        visNetwork::visEdges(color = "#A9A9A9", width = 2) |>
        visNetwork::visIgraphLayout(layout)
}



#' @return A plot
#'
#' @export
#' @import plotly
#' @rdname pep_prot_CC
#'
plotCCJitter <- function(
        df) {

    # i_tooltip <- which(startsWith(colnames(df), "tooltip"))
    txt_tooltip <- NULL

    # if (length(i_tooltip) == 0){
    #  warning("There is no tooltip in the object.")
    # } else {
    # for (i in i_tooltip) {
    #   txt_tooltip <- paste(txt_tooltip, "<b>",
    #                        gsub("tooltip_", "", colnames(df)[i], fixed = TRUE),
    #                        " </b>: {point.", colnames(df)[i], "} <br> ",
    #                        sep = "")
    # }
    # }

    p <- plot_ly(
        data = df,
        x = ~x,
        y = ~y,
        type = "scatter",
        mode = "markers",
        text = txt_tooltip,
        hoverinfo = "text",
        showlegend = FALSE
    ) |>
        plotly::layout(
            xaxis = list(title = "Nb of peptides ic CC"),
            yaxis = list(title = "Nb of proteins ic CC"),
            hovermode = "closest",
            margin = list(b = 60) 
        )
    
    return(p)
}





#' @title Connected Components infos
#' @param cc A list of connected component
#' @return A `list` of three items:
#' * `One_One`: the number of cc composed of one protein and one peptide
#' * `One_Multi`: the number of cc composed of one protein and several peptides
#' * `Multi_Multi`: the number of cc composed of several proteins and
#' several (shared) peptides.
#'
#' @examples
#' if (interactive()) {
#'     data(sub_R25)
#'     GetCCInfos(get_cc(sub_R25[[1]]))
#' }
#'
#' @export
#' @rdname pep_prot_CC
#'
GetCCInfos <- function(cc) {
    # stopifnot(inherits(cc, "list"))
    cc.infos <- list(
        One_One = list(),
        One_Multi = list(),
        Multi_Multi = list()
    )


    ll.prot <- lapply(cc, function(x) {
        ncol(x)
    })
    ll.pept <- lapply(cc, function(x) {
        nrow(x)
    })
    ll.prot.one2one <- intersect(which(ll.prot == 1), which(ll.pept == 1))
    ll.prot.one2multi <- intersect(which(ll.prot == 1), which(ll.pept > 1))
    ll.prot.multi2any <- which(ll.prot > 1)

    cc.infos[["One_One"]] <- cc[ll.prot.one2one]
    cc.infos[["One_Multi"]] <- cc[ll.prot.one2multi]
    cc.infos[["Multi_Multi"]] <- cc[ll.prot.multi2any]

    cc.infos
}
