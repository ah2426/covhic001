#' visualizes conditional independence network
#' 
#' NOTE
#' 1. highlight cDC-CCL13 axis
#' 2. edge color by marginal spearman correlation
#'
#' Name: Ao H
#' Date: 2025-09-10

rm(list = ls())
library(tidyverse)
library(igraph)
library(gplots)
library(ggraph)

# load data
MATR_CORR = here("data", "baseline_spearman_correlation_object.RDS")
HUGE_FIT  = here("outs", "network", "huge_model.RDS")

rho_corr_object = readRDS(MATR_CORR)
data = rho_corr_object$r # correlation matrix
est = readRDS(HUGE_FIT)

# save dir
SAVE_DIR = here("outs", "network")
dir.create(SAVE_DIR, showWarnings = FALSE)

# helper function
source(here("R/src", "plt_correlation_heatmap.R"))


# VISUAL CORRELATION HEATMAP -------------------------------------------------
cat('\n\n',str_pad(' VISUAL CORRELATION HEATMAP ',80,'both','-'),'\n\n',sep='')

p = plt_correlation_heatmap(
    corr = rho_corr_object$r,
    pval = rho_corr_object$P,
    legend.title = "Spearman", title = ""
)

pdf(file.path(SAVE_DIR, "efig7d-spearman_correlation_heatmap.pdf"), 10, 10)
print(p)
dev.off()


# VISUAL NETWORK -------------------------------------------------------------
cat('\n\n',str_pad(' VISUAL NETWORK ',80,'both','-'),'\n\n',sep='')

# target nodes - CCL13-cDC axis
TARGET_NODES = c("MCP-4",
                 "CD1c+_DC_freq._of_total",
                 "CD1c-_DC_freq._of_total")

# (opt) rename some variables
feature_rename = list(
    `MCP-4` = "CCL13", 
    `MDC` = "CCL22",
    `CD1c+_DC_freq._of_total` = "CD1c+ DCs",
    `CD1c-_DC_freq._of_total` = "CD1c- DCs",
    `Ratio_CD57−/CD57+` = "Ratio CD56dim CD57-/CD57+ NK",
    `Sum_of_IL-2_NSP_pool_response` = "Sum of RTC pools IL-2 response",
    `CD56bright_Freq._of_Parent` = "CD56bright NK"
)
colnames(data) = dplyr::recode(colnames(data), !!!feature_rename)
TARGET_NODES = dplyr::recode(TARGET_NODES, !!!feature_rename)
rownames(data) = colnames(data)

# ADJ MATRIX OF SELECTED LAMBDA -------------

opt_index = est$opt.index
opt_lambda = est$opt.lambda
cat("\n--> Optimal lambda:", opt_lambda, "\n")
cat("--> Optimal lambda index:", opt_index, "\n")

adj_matrix = est$path[[opt_index]]
colnames(adj_matrix) = colnames(data)
rownames(adj_matrix) = colnames(data)


# GRAPH STRUCTURE --------------------

graph = graph_from_adjacency_matrix(
    adj_matrix, 
    mode = "undirected"
)


# TARGET & NEIGHBORS --------------------

# target nodes and their 1st-order neighbors
target_index = which(V(graph)$name %in% TARGET_NODES)
neighbor_index = unique(unlist(ego(graph, order = 1, nodes = TARGET_NODES)))

# emphasize target nodes & 1st-order neighbors
V(graph)$alpha = ifelse(seq_along(V(graph)) %in% neighbor_index, 1, 0)
V(graph)$show_label = ifelse(seq_along(V(graph)) %in% neighbor_index, TRUE, FALSE)
E(graph)$alpha = apply(ends(graph, E(graph)), 1, function(x) {
    if (all(x %in% V(graph)$name[neighbor_index])) 1 else 0
})


# NODE TYPE --------------------

# add node type
V(graph)$node_type = dplyr::case_when(
    names(V(graph)) == "Sum of RTC pools IL-2 response" ~ "T cell response",
    names(V(graph)) %in% c("CD1c- DCs", "CD1c+ DCs", 
        "CD56bright NK", "Ratio CD56dim CD57-/CD57+ NK") ~ "Cell frequency",
    TRUE                                                 ~ "Cytokine"
)

# add corresponding colors to node_type
my_colors = c(
    "Cytokine" = "gold",
    "Cell frequency" = "tomato",
    "T cell response" = "gray40"
)

# annotate edge with correlation
for (edge in E(graph)) {
    
    nodes = ends(graph, edge)
    node1 = nodes[1]
    node2 = nodes[2]
    
    # get spearman correlation from "data" matrix
    corr = data[node1, node2]
    
    if (all(nodes %in% V(graph)$name[neighbor_index])) {
        E(graph)[edge]$rho_corr = corr
    } else {
        E(graph)[edge]$rho_corr = NA_real_
    }
}


# ggplot version of graph using ggraph --------------------

layout = layout_with_kk(graph)
V(graph)$x = layout[,1]
V(graph)$y = layout[,2]

df_highlight <- tibble(
    name = V(graph)$name,
    x = V(graph)$x, y = V(graph)$y,
    is_target = V(graph)$name %in% TARGET_NODES
)

p = ggraph(graph, layout = "manual", x = V(graph)$x, y = V(graph)$y) +
    # convex hull shading around target group
    ggforce::geom_mark_hull(
        data = df_highlight,
        aes(filter = is_target, x = x, y= y),
        concavity = 20, expand = unit(5, "mm"),
        alpha = 0.2, fill = "grey90", colour = "grey70"
    ) +
    
    # edges
    geom_edge_link(
        aes(color = rho_corr, alpha = alpha), # , width = abs(rho_corr)
        width = 2,
    ) +
    scale_edge_color_gradient2(
        low = "#1F78B4", mid = "grey100", high="black",
        midpoint = 0,
        name = "Spearman",
        na.value = NA
    ) + 
    scale_edge_alpha(range = c(0, 1), guide = "none") +

    # nodes
    geom_node_point(
        aes(fill = node_type, alpha = alpha), 
        size = 12, color = "black", shape = 21
    ) +
    geom_node_text(
        aes(label = ifelse(show_label, name, "")), 
        size = 5, repel = FALSE
    ) +
    scale_fill_manual(
        values = my_colors, 
        name = "Node type"
    ) +
    scale_alpha_identity() +  # use alpha values as-is (don’t add legend)
    theme_void(base_size=16) +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)
    )

pdf(file.path(SAVE_DIR, "fig6c-network_ccl13_cdc_axis_ggraph.pdf"), width=10, height=8)
print(p)
dev.off()


# ANALYSIS FINISHED ----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
