#' visualizes conditional independence network
#' 
#' NOTE
#' 1. highlight cDC-CCL13 axis
#' 2. edge weighted by partial correlation (or Spearman correlation)
#' 3. partial correlation obtained by transforming the estimated precision matrix
#'
#' Name: Ao H
#' Date: 2025-06-25
#' Envs: Seurat5

re_init()
library(tidyverse)
library(igraph)
library(gplots)
library(ggraph)

DATA_IDX = "lst_rho_corr-all-use_copula_FALSE"
PD_PROJ = FALSE
TARGET_NODES = c("MCP-4",
                 "CD1c+_DC_freq._of_total",
                 "CD1c-_DC_freq._of_total") # "IL-10", "IL-18"

# DIR CONFIG --------------------
WDIR = "~/projects/muco"
PROJ = file.path(WDIR,"scripts/network_analysis")

SAVE_DIR = file.path(PROJ, "outs/huge_fit", DATA_IDX,
                     paste0("pd_proj_", PD_PROJ))
mkdir(SAVE_DIR)

SCRIPT_DIR = file.path(PROJ, "src")
source_all(SCRIPT_DIR)

# input file paths
EST = file.path(SAVE_DIR, "model.RDS")
DATA= file.path(SAVE_DIR, "rho_corr_matrix.RDS")

# load data
est = readRDS(EST)
data = readRDS(DATA)

# (opt) rename some the variables
feature_rename = list(
    `MCP-4` = "CCL13", `MDC` = "CCL22",
    `CD1c+_DC_freq._of_total` = "CD1c+ DCs %",
    `CD1c-_DC_freq._of_total` = "CD1c- DCs %",
    `Ratio_CD57−/CD57+` = "CD56dim:CD57-/CD57+ NK",
    `Sum_of_IL-2_NSP_pool_response` = "Sum of IL-2 NSP pool response",
    `CD56bright_Freq._of_Parent` = "CD56bright NK %"
)
colnames(data) = dplyr::recode(colnames(data), !!!feature_rename)
TARGET_NODES = dplyr::recode(TARGET_NODES, !!!feature_rename)
rownames(data) = colnames(data)

# ADJ MATRIX OF SELECTED LAMBDA -------------

print(est$lambda)
opt_index = est$opt.index
opt_lambda = est$opt.lambda
cat("\n--> Optimal lambda:", opt_lambda, "\n")
cat("--> Optimal lambda index:", opt_index, "\n")

adj_matrix = est$path[[opt_index]]
colnames(adj_matrix) = colnames(data)
rownames(adj_matrix) = colnames(data)

    
# GRAPH STRUCTURE --------------------

# igraph
graph = graph_from_adjacency_matrix(
    adj_matrix, 
    mode = "undirected"
)

# Target node and its neighbor_indexs
target_index = which(V(graph)$name %in% TARGET_NODES)
neighbor_index = unique(unlist(ego(graph, order = 1, nodes = TARGET_NODES)))
print(neighbor_index)

# add node type
V(graph)$node_type = dplyr::case_when(
    names(V(graph)) == "Sum of IL-2 NSP pool response" ~ "T cell response",
    names(V(graph)) %in% c("CD1c- DCs %", "CD1c+ DCs %", 
                           "CD56bright NK %", "Ratio_CD57-/CD57+") ~ "Cell frequency",
    TRUE ~ "Cytokine"
)

# add corresponding colors to node_type
my_colors = c(
    "Cytokine"        = "gold",
    "Cell frequency"  = "tomato",
    "T cell response" = "gray40"
)


# default graph display setting ----------------------

V(graph)$color = NA                    # no colors
V(graph)$size  = 2
V(graph)$label.cex = 1.5
V(graph)$label.color = "black"
V(graph)$frame.color = NA              # no node.frame.colors
V(graph)$frame.width = 1
E(graph)$color = NA
E(graph)$width = 0

# only show labels for targeted nodes & neighbors
vertex_labels = rep(NA, vcount(graph)) # no labels
vertex_labels[neighbor_index] = V(graph)$name[neighbor_index]
print(vertex_labels)

# highlight target node & 1st order neighbor
V(graph)$color[neighbor_index] = my_colors[V(graph)$node_type[neighbor_index]]
V(graph)$frame.color[neighbor_index] = "black"
V(graph)$frame.width[neighbor_index] = 1
V(graph)$size[neighbor_index]  = 8


# set edge color & transparency based on correlation --------------------

print(neighbor_index) # include target nodes
edges_pairwise_neighbor = combn(neighbor_index, 2, simplify = FALSE)
for (pair in edges_pairwise_neighbor) {
    
    # if edge exist between a pair
    if (are_adjacent(graph, pair[1], pair[2])) {
        edge = E(graph, P = pair)
        node1 = ends(graph, edge)[1]
        node2 = ends(graph, edge)[2]
        
        # get spearman correlation from "data" matrix
        # partial correlation only applicable using glasso
        corr = data[node1, node2]
        
        # Set color based on sign
        E(graph)[edge]$color = if (corr > 0) "black" else "#1F78B4"
        
        # Set transparency based on corr level
        E(graph)[edge]$corr = corr
        
        # highlight edge width based on correlation
        E(graph)[edge]$width = 10 * abs(corr)
    }
}


# # set edge transparency based on spearman correlation --------------------
# # set weight to 0 if NA
# E(graph)$weight[is.na(E(graph)$weight)] = 0
# 
# edge_alpha = scales::rescale(
#     abs(E(graph)$weight), 
#     to = c(0.3, 1) # 0.1 = faint, 1 = opaque
# )

# mark nodes --------------------

mark_groups = list(target_index) # the CCL13-cDC axis


# graph layout --------------------

#layout = layout_with_fr(graph)
layout = layout_with_kk(graph)
font_family = "Helvetica"

df_map = tibble(width = E(graph)$width, corr = abs(E(graph)$corr)) |>
    drop_na() |> arrange(desc(corr)) |>
    print()

corr_min = min(df_map$corr)
corr_max = max(df_map$corr)

# 5 evenly spaced correlation values
legend_corr = seq(corr_min, corr_max, length.out = 4)

# interpolate widths for those values
legend_width = approx(x = df_map$corr, y = df_map$width, xout = legend_corr)$y

pdf(file.path(SAVE_DIR, "network_ccl13_cdc_axis_edge_width.pdf"), width=12, height=12)

plot(
    graph, layout = layout,
    vertex.label = vertex_labels,
    size = 4,
    vertex.size = V(graph)$size, 
    vertex.label.cex = V(graph)$label.cex, 
    vertex.label.color = V(graph)$label.color,
    vertex.frame.color = V(graph)$frame.color, 
    vertex.frame.width = V(graph)$frame.width,
    vertex.label.family = font_family,
    vertex.label.dist = 1,     # distance from center of node
    vertex.label.degree = 3*pi/2,
    vertex.color = V(graph)$color, 
    edge.width = E(graph)$width, 
    edge.color = E(graph)$color,
    # edge.color = rgb(1,1,1, alpha=edge_alpha), # edge_alpha based on correlation
    # edge.curved=.1, # add curve
    mark.groups = mark_groups,
    mark.col=c("#F2F2F2") , mark.border = "grey60"
)

# legend of node class
legend(
    "topleft",
    legend = names(my_colors),
    pch    = 21,
    col    = "white",
    pt.bg  = my_colors,
    pt.cex = 4,
    cex    = 1.5,
    bty    = "n", # no box
    title  = "Node type"
)

# legend of spearman correlation values
legend(
    "bottomleft",
    legend = sprintf("%.2f", legend_corr),
    lwd    = legend_width,
    cex    = 1.5,
    bty    = "n",
    title  = "Spearman"
)

dev.off()


# ANALYSIS FINISHED ----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
