#' wrangles/normalizes baseline feature data for elastic-net & network analysis
#' 
#' INPUT
#' 1. nasal soluble mediators
#' 2. cell frequencies
#' 3. T cell response (fluorospot)
#' 
#' Name: Ao H
#' Date: 2025-09-10

rm(list = ls())
library(tidyverse)

# load data
df_cfreq = read_csv(here("data/data_original", "cell_frequencies.csv")) |>
    mutate(Participant_ID = as.character(Participant_ID))
df_fluorospot = read_csv(here("data/data_original", "fluorospot.csv")) |>
    mutate(Participant_ID = as.character(Participant_ID))
df_nasal = read_csv(here("data/data_original", "nasal_soluble_mediators.csv")) |>
    mutate(Participant_ID = as.character(Participant_ID))


# prep data for elastic-net training --------------------

# log10-stablize & scale
df_nasal_log = df_nasal |> mutate(across(where(is.numeric), ~ log10(.)))
df_nasal_norm = df_nasal_log |> mutate(across(where(is.numeric), ~ as.numeric(scale(.))))

write_csv(df_nasal_norm, here("data", "nasal_soluble_mediators_normalized.csv"))


# prep data for network inference --------------------

# feature rename
colnames(df_cfreq) = colnames(df_cfreq) |> str_replace_all(" ", "_")
colnames(df_fluorospot) = colnames(df_fluorospot) |> str_replace_all(" ", "_")

# merge nasal soluble proteins, cell frequencies, T cell response
df_merge = df_nasal_log |>                                             
    full_join(df_cfreq, by=c("Participant_ID", "Group")) |>
    full_join(df_fluorospot, by=c("Participant_ID", "Group"))

# keep only the ratio CD56dim:CD57-/CD57+
df_merge = df_merge |>
    select(-c(`CD56dim_CD57-_Freq._of_Parent`, `CD56dim_CD57+_Freq._of_Parent`))

write_csv(df_merge, file.path(here("data", "baseline_features_merge.csv")))

# compute pairwise spearman correlation matrix 
# as input for huge to infer condi independence network
# ** auto-remove NA entries

X = df_merge |> select(-c(Group, Participant_ID)) |> as.matrix()
rho_corr_object = Hmisc::rcorr(X, type="spearman")
print(names(rho_corr_object))

saveRDS(rho_corr_object,
        file.path(here("data", "baseline_spearman_correlation_object.RDS")))


# ANALYSIS FINISHED -----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
