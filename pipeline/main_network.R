#' estimates conditional independence network using huge
#' 
#' NOTE
#' 1. project to the nearest PD matrix or not
#' 
#' Name: Ao Huang
#' Date: 2025-03-18 (modified 2025-08-28)
#' Envs: Seurat5

re_init()
library(tidyverse)
library(huge)
set.seed(72)

DATA_IDX = "lst_rho_corr-all-use_copula_FALSE"
PD_PROJ = TRUE
TARGET_NODE = "MCP-4"

# DIR CONFIG --------------------
WDIR = "~/projects/muco"
PROJ = file.path(WDIR,"scripts/network_analysis")

SAVE_DIR = file.path(PROJ, "outs/huge_fit", DATA_IDX,
                     paste0("pd_proj_", PD_PROJ))
mkdir(SAVE_DIR)

SCRIPT_DIR = file.path(PROJ, "src")
source_all(SCRIPT_DIR)

# input file paths
INPUT_MATR = file.path(PROJ, "outs/rho_matrix", paste0(DATA_IDX, ".RDS"))

# # load data
data = readRDS(INPUT_MATR)
data = data$corr


# PDS_PROJECTION ------------------------------------------------------------
cat('\n\n',str_pad(' PDS_PROJECTION ',80,'both','-'),'\n\n',sep='')

if (PD_PROJ) {
    message("... projecting to the nearest P.D. matrix")
    near_pd_object = Matrix::nearPD(data, corr=TRUE)
    data = near_pd_object$mat |> as.matrix()
} else {
    message("... FALSE: projection nearest P.D. matrix")
}


# HUGE FIT -------------------------------------------------------------------
cat('\n\n',str_pad(' HUGE FIT ',80,'both','-'),'\n\n',sep='')

est = huge(
    data, 
    method    = "mb", # graphical lasso
    nlambda   = 100,
    sym       = "or"
)

# select the opt. model -------------------

# OPTION 1
# find elbow point of the sparsity-lambda curve
# ** can also try to implement a CV framework

df_fit = tibble(sparsity = est$sparsity, lambda = est$lambda)
opt_lambda_idx = pathviewr::find_curve_elbow(df_fit)
opt_lambda = df_fit$lambda[opt_lambda_idx]
opt_lambda_idx = which(est$lambda == opt_lambda)

# OPTION 2 - use ebic
# Model selection is not available when using the covariance matrix as input
# fit_opt = huge.select(est, criterion = "ebic", ebic.gamma = 0.5)

# OPTION 3 - compute ebic for each huge fit
# only applicable if method = "glasso"
# ebic_path = ebic_glasso(
#     huge_fit = est,
#     n_sample = 34, # 16 uninfected 18 infected
#     gamma = 1e-1,
#     near_zero = 1e-12
# )
# est$ebic_path = ebic_path
# plot(est$lambda, est$ebic_path, type="b", xlab="log(lambda)", ylab="EBIC")
# opt_lambda_idx = which.min(ebic_path)

opt_lambda = est$lambda[opt_lambda_idx]
cat("--> Optimal lambda:", opt_lambda, "\n")
cat("--> Optimal lambda index:", opt_lambda_idx, "\n")
est$opt.lambda = opt_lambda
est$opt.index = opt_lambda_idx


# report primary nodes to the target node --------------

cat("\n--> Primary neighbors of", TARGET_NODE, ":\n")

adj_matrix = est$path[[opt_lambda_idx]]
colnames(adj_matrix) <- rownames(adj_matrix) <- colnames(data)
primary_nodes = rownames(adj_matrix)[adj_matrix[TARGET_NODE,] != 0]
print(primary_nodes)

# save huge fit
saveRDS(est, file.path(SAVE_DIR, "model.RDS"))
saveRDS(data, file.path(SAVE_DIR, "rho_corr_matrix.RDS"))


# HUGE FIT QC ----------------------------------------------------------------
cat('\n\n',str_pad(' HUGE FIT QC ',80,'both','-'),'\n\n',sep='')

cat('--> lambda path', '\n')
cat('--> opt. lambda', '\n')
pdf(file.path(SAVE_DIR, "lambda_path.pdf"), width = 6, height = 4)
    print(plt_opt_lambda(est, opt_index = opt_lambda_idx))
    plot(est)
dev.off()

cat('\n--> feature path', '\n')
pdf(file.path(SAVE_DIR, "feature_path.pdf"), width = 8, height = 6)
    p = plt_feature_path(
        est = est,
        data = data,
        target = TARGET_NODE,
        opt_index = opt_lambda_idx
    )
    print(p)
dev.off()


# ANALYSIS FINISHED -----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
