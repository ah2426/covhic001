#' main script to run elastic-net regression (family: binomial)
#'
#' Name: Ao H
#' Date: 2025-09-10

rm(list = ls())
random_seed = 42
set.seed(random_seed) # random seed for reproducibility

library(tidyverse)
library(here)
library(eNetXplorer)
library(future)
library(furrr)  # parallel

# Set up parallel processing
N_CORES = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1L))  
future::plan(future::multicore, workers=N_CORES)

# CONSOLE INPUT -------------------

# data file
DATA = here("data", "nasal_soluble_mediators_normalized.csv")

# save dir
SAVE_DIR = here("outs/elastic-net")
dir.create(SAVE_DIR, recursive=TRUE, showWarnings=FALSE)
setwd(SAVE_DIR)

# helper functions
SCRIPT = file.path(PROJ,'src')
source_all(SCRIPT)

# TRAINING CONFIG --------------------
# ** configuration before run

FAMILY = "binomial"
ALPHAS = seq(0,1,by=0.05) # 0:ridge, 1:lasso
N_RUN  = 100   # number of cross-validation runs
N_PERM = 100   # number of permutation to generate null distribution
N_FOLD = 10    # CV fold
METRIC = "acc" # select model with the best accuracy

cat('\n# CONSOLE REPORT ----------------\n')
cat('DATA_FILE:', DATA, '\n')
cat('SAVE_DIR:', SAVE_DIR, '\n')
cat('FAMILY:', FAMILY, '\n')
cat('ALPHAS:', ALPHAS, '\n')
cat('N_PERM:', N_PERM, '\n')
cat('N_RUN:', N_RUN, '\n')
cat('N_FOLD:', N_FOLD, '\n')
cat('METRIC:', METRIC, '\n')
cat('N_CORES (PARALLEL):', N_CORES, '\n')


# LOAD & PREPROC DATA --------------------------------------------------------
cat('\n\n',str_pad(' LOAD & PREPROC DATA ',80,'both','-'),'\n\n',sep='')

# load data
df = read_csv(DATA, show_col_types = FALSE) |> 
    as.data.frame() |> column_to_rownames("Participant_ID") |>
    mutate(Group = ifelse(Group == "Uninfected", 1, 0)) # binarize
head(df)

# for prediction
y = df |> select(Group) |> as.matrix() # 0:Infected, 1:Uninfected
X = df |> select(-Group) |> as.matrix()

cat('--> data dim:', '\n')
print(dim(df))


# TRAIN ELASTIC NET -----------------------------------------------------------
cat('\n\n',str_pad(' TRAIN ELASTIC NET ',80,'both','-'),'\n\n',sep='')

toc = Sys.time()

# # ** parallel run
eNet = ALPHAS %>% furrr::future_map(
    ~ eNetXplorer(
        x=X, y=y,
        family=FAMILY,
        scaled = FALSE, # scaled is done/or_not to data input already
        binom_method = METRIC,
        n_run=N_RUN, 
        n_perm_null=N_PERM, 
        n_fold=N_FOLD, seed=42,
        save_obj=T, dest_obj=paste0("eNet_a",.x,".Robj"),
        alpha=.x
    ),
    .options = furrr_options(seed=TRUE)
)

# save family of enet models
mergeObj( paste0("eNet_a", ALPHAS, ".Robj") )
load("./eNet_merged.Robj")
cat('\n--> eNet model SAVED \n')

tic = Sys.time()

cat('\n--> Time Spent:\n')
print(tic-toc)

cat('\n--> Model Summary:\n')
print(summary(eNet))


# ANALYSIS FINISHED -----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
