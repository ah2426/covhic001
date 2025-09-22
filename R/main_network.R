#' main script to estimate conditional independence network using huge
#' 
#' Name: Ao H
#' Date: 2025-09-10

rm(list = ls())
set.seed(72)
library(tidyverse)
library(here)
library(huge)

# load data
# correlation matrix as huge input
MATR_CORR = here("data", "baseline_spearman_correlation_object.RDS")
rho_corr_object = readRDS(MATR_CORR)
data = rho_corr_object$r

# save dir
SAVE_DIR = here("outs", "network")
dir.create(SAVE_DIR, showWarnings = FALSE)

# helper function
source(here("R/src", "plt_opt_lambda.R"))


# HUGE FIT -------------------------------------------------------------------
cat('\n\n',str_pad(' HUGE FIT ',80,'both','-'),'\n\n',sep='')

est = huge(
    data, 
    method  = "mb",
    nlambda = 100,
    sym     = "or"
)

# select the opt. model -------------------
# find elbow point of the sparsity-lambda curve

df_fit = tibble(sparsity = est$sparsity, lambda = est$lambda)
opt_lambda_idx = pathviewr::find_curve_elbow(df_fit)
opt_lambda = df_fit$lambda[opt_lambda_idx]
opt_lambda_idx = which(est$lambda == opt_lambda)

opt_lambda = est$lambda[opt_lambda_idx]
cat("--> Optimal lambda:", opt_lambda, "\n")
cat("--> Optimal lambda index:", opt_lambda_idx, "\n")
est$opt.lambda = opt_lambda
est$opt.index = opt_lambda_idx

# lambda path plots
cat('--> lambda path', '\n')
cat('--> opt. lambda', '\n')
plt_opt_lambda(est, opt_index = opt_lambda_idx)
plot(est)

# save huge fit
saveRDS(est, file.path(SAVE_DIR, "huge_model.RDS"))


# ANALYSIS FINISHED -----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
