#' visualizes elastic-net model training output
#'
#' Name: Ao H
#' Date: 2025-09-10

rm(list = ls())
library(tidyverse)
library(here)
library(eNetXplorer)
library(yardstick)

# DIR CONFIG --------------------

ENET_MODEL = here("outs/elastic-net", "eNet_merged.Robj")
META = here("data", "meta_subject.csv")

# save dir
SAVE_DIR = here("outs/elastic-net")

# helper functions
SCRIPT_DIR = here("R", "src")
source(file.path(SCRIPT_DIR, "plt_error_bar.R"))
source(file.path(SCRIPT_DIR, "plt_roc_curve.R"))

# load data
load(ENET_MODEL) # eNet
meta_subject = read_csv(META, show_col_types = TRUE) |>
    drop_na() |> mutate(Participant_ID = as.character(Participant_ID)) |>
    rename(true.id = Group, true.id.2 = Group.2) |>
    print()


# MODEL PERFORMANCE ----------------------------------------------------------
cat('\n\n',str_pad(' MODEL PERFORMANCE ',80,'both','-'),'\n\n',sep='')

# best lambda (best pvalue)
alpha_optimal = which.max(-log10(eNet$QF_model_vs_null_pval))

# prediction accuracy across alpha runs
pdf(file.path(SAVE_DIR, "efig7a-train_summary.pdf"), 4.5, 4.5)
plot.eNetXplorer(eNet, plot.type="summary")
dev.off()

df_fit = .prep_caterpillar(
    x=eNet, 
    alpha.index=alpha_optimal, 
    n_top_features=5, 
    stat="freq"
)
top_features = df_fit$feature

# feature frequency permutation significance
p1 = plt_error_bar(
    x = eNet,
    alpha.index = alpha_optimal,
    stat = "freq",
    n_top_features = 5,
    features = top_features
) + theme(legend.position = "none") +
    labs(y = "frequency")

# coefficient permutation significance
p2 = plt_error_bar(
    x = eNet,
    alpha.index = alpha_optimal,
    stat = "coef",
    n_top_features = 5,
    features = top_features
) + theme(legend.position = "none") + 
    theme(axis.text.y = element_blank(),
          plot.title = element_blank()) +
    labs(y = "coefficient")

library(patchwork)
p = p1 + p2
pdf(file.path(SAVE_DIR, "fig6b-enet_top_features.pdf"), 7, 4)
print(p)
dev.off()


# PREDICTION PROB ------------------------------------------------------------
cat('\n\n',str_pad(' PREDICTION PROB ',80,'both','-'),'\n\n',sep='')

# probability across repeated CV runs that is assigned the "Uninfected" label
# 0:Infected, 1:Uninfected
df_prob = eNet$predicted_values[[alpha_optimal]] |>
    as.data.frame() |> select(`1`) |> rename(probability_uninfected=`1`) |>
    rownames_to_column("Participant_ID") |> as_tibble() |>
    print()

# include group information
df_prob = df_prob |>
    inner_join(meta_subject, by="Participant_ID") |>
    mutate(.POS_probability = probability_uninfected) |>
    mutate(true.id = factor(true.id, levels=c("Infected", "Uninfected"))) |>
    mutate(true.id.2 = recode(true.id.2, "Sustained" = "Infected")) |>
    mutate(true.id.2 = factor(true.id.2, levels = c("Infected", "Transient", "Abortive"))) |>
    print()

# prob distribution --------------------

p = ggplot(df_prob, aes(x=true.id, y=probability_uninfected)) +
    geom_boxplot(
        aes(color = true.id),
        outlier.shape = NA, 
        width=.5,
        show.legend = FALSE
    ) + 
    scale_color_manual(values=c(Uninfected = "blue", Infected = "red")) +
    
    ggnewscale::new_scale_color() + # reset color scale 
    
    ggbeeswarm::geom_beeswarm(
        aes(color=true.id.2, shape=true.id.2), 
        size=3, cex=4, stroke=1
    ) +
    scale_shape_manual(values=c("Infected"=1, "Abortive"=0, "Transient"=2)) +
    scale_color_manual(values=c("Infected"="red", "Abortive"="blue", "Transient"="darkgrey")) +
    
    ggpubr::stat_compare_means(
        # aes(group = true.id)
        comparisons = list(c("Uninfected","Infected"))
    ) + 
    theme_classic(base_size=16) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    #theme(axis.text.x = element_blank()) +
    labs(x="", y="Probability of Uninfected", )
pdf(file.path(SAVE_DIR, "efig7b-predict_probability_across_runs.pdf"), 5, 5)
print(p)
dev.off()


# AUROC curve --------------------

df_prob = df_prob |>
    mutate(true.id = factor(true.id, levels=c("Uninfected", "Infected")))

CLASS = "Uninfected"
auroc = roc_auc(df_prob, truth=true.id,.POS_probability)
df_roc = roc_curve(df_prob, truth=true.id, .POS_probability)

# auroc curve
auroc_score = round(auroc$.estimate, 3)
p = plt_roc_curve(df_roc, auroc_value = auroc_score)
pdf(file.path(SAVE_DIR, "fig6a-auroc_across_runs.pdf"), 4, 4)
print(p)
dev.off()


# plot distribuion of IL-10 --------------------

DF_RAW = here("data", "nasal_soluble_mediators_raw.csv")
df_protein = read_csv(DF_RAW, show_col_types = FALSE) |>
    mutate(Participant_ID = as.character(Participant_ID)) |>
    inner_join(meta_subject, by="Participant_ID") |>
    mutate(true.id = factor(true.id, levels=c("Infected", "Uninfected"))) |>
    mutate(true.id.2 = recode(true.id.2, "Sustained" = "Infected")) |>
    mutate(true.id.2 = factor(true.id.2, levels = c("Infected", "Transient", "Abortive"))) |>
    print()

p = ggplot(df_protein, aes(x=true.id, y=`IL-10`)) +
    geom_boxplot(
        aes(color = true.id),
        outlier.shape = NA, 
        width=.5,
        show.legend = FALSE
    ) + 
    scale_color_manual(values=c(Uninfected = "blue", Infected = "red")) +
    
    ggnewscale::new_scale_color() + # reset color scale 
    
    ggbeeswarm::geom_beeswarm(
        aes(color=true.id.2, shape=true.id.2), 
        size=3, cex=4, stroke=1
    ) +
    scale_shape_manual(values=c("Infected"=1, "Abortive"=0, "Transient"=2)) +
    scale_color_manual(values=c("Infected"="red", "Abortive"="blue", "Transient"="darkgrey")) +
    
    ggpubr::stat_compare_means(
        # aes(group = true.id)
        comparisons = list(c("Uninfected","Infected"))
    ) + 
    
    scale_y_log10() +
    
    theme_classic(base_size=16) + 
    theme(axis.text.x = element_text(angle=45, hjust=1)) + 
    labs(x="", y="IL-10 (pg/ml)")

pdf(file.path(SAVE_DIR, "efig7c-log_il10_distribution.pdf"), 5, 5)
print(p)
dev.off()


# ANALYSIS FINISHED -----------------------------------------------------------
cat('\n\n',str_pad(' ANALYSIS FINISHED ',80,'both','-'),'\n\n',sep='')
