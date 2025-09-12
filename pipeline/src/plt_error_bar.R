# elastic-net related plot functions

plt_error_bar = function (x, # eNet object
                          alpha.index, 
                          n_top_features=10, 
                          features = NULL,
                          title = "", 
                          stat="freq") {
    
    # preps plot dataframe
    plt_df = .prep_caterpillar(
        x = x,
        alpha.index = alpha.index,
        n_top_features = n_top_features,
        features, 
        stat = stat
    ) |>
        pivot_longer(
            c(feature_mean, null_feature_mean),
            names_to = "Type", values_to = "mean"
        ) |>
        mutate(
            lower = ifelse(Type == "feature_mean", lower_bound, lower_bound_null),
            upper = ifelse(Type == "feature_mean", upper_bound, upper_bound_null),
            Type  = ifelse(Type == "feature_mean", "Model", "Null")
        ) |>
        mutate(Type = factor(Type, levels = c("Null", "Model"))) |>
        mutate(pval_star = ifelse(Type == "Model", pval_star, ""))
    
    # show significance
    df_pval = plt_df |> filter(Type == "Model", pval < .05) |>
        mutate(pos = ifelse(mean >= 0, mean + 0.6, mean - 0.6))
    
    p = ggplot(plt_df, aes(x = reorder(feature, desc(pval)), y = mean, fill = Type)) +
        geom_col(width = 0.6, position = position_dodge(width=0.6)) +
        geom_errorbar(aes(ymin = lower, ymax = upper),
                      width = 0.3, position = position_dodge(width=0.6)) +
        
        # add sig level
        geom_text(
            data = df_pval,
            aes(label = pval_star, y = pos),
            vjust = 0.5
        ) +
        
        coord_flip() +
        labs(x = "", y = "", fill = "Type") +
        scale_fill_manual(values = c(Model="#D95F02", Null="grey")) +
        ggtitle(title) +
        theme_classic(base_size=16)
    
    # for coef - symmetric x-axis
    if (stat == "coef") {
        max_abs = max(abs(c(plt_df$upper, plt_df$lower)), na.rm = TRUE) + 0.5
        p = p + ylim(-max_abs, max_abs) +
            geom_hline(yintercept = 0, color = "black") + 
            theme(axis.line.y = element_blank(),
                  axis.ticks.y = element_blank())
    }
    
    return(p)
}


.prep_caterpillar = function(x, alpha.index, n_top_features, stat, features = NULL) {
    
    if (stat == "coef") {
        df = data.frame(
            feature_mean = x$feature_coef_wmean[, alpha.index],
            feature_sd   = x$feature_coef_wsd[, alpha.index],
            null_feature_mean = x$null_feature_coef_wmean[, alpha.index],
            null_feature_sd   = x$null_feature_coef_wsd[, alpha.index],
            pval = x$feature_coef_model_vs_null_pval[, alpha.index]
        )
    } else if (stat == "freq") {
        df = data.frame(
            feature_mean = x$feature_freq_mean[, alpha.index],
            feature_sd   = x$feature_freq_sd[, alpha.index],
            null_feature_mean = x$null_feature_freq_mean[, alpha.index],
            null_feature_sd   = x$null_feature_freq_sd[, alpha.index],
            pval = x$feature_freq_model_vs_null_pval[, alpha.index]
        )
    }
    
    df = df |> rownames_to_column("feature") |>
        filter(!if_any(everything(), is.nan)) |> # remove NaN rows
        arrange(pval)
    
    # restrict to selected features
    if (!is.null(features)) {
        df = df |> filter(feature %in% features)
    } else {
        df = df |> head(n_top_features)
    }
    
    # add sig level to feature names
    df$pval_star = ""
    df$pval_star[df$pval < .05]  = "*"
    df$pval_star[df$pval < .01]  = "**"
    df$pval_star[df$pval < .001] = "***"
    #df$feature = paste(df$feature, df$pval_star)
    
    # freq-specific lower bound
    if (stat == "freq") {
        df = df |> mutate(
            lower_bound = pmax(0, feature_mean - feature_sd),
            upper_bound = pmin(1, feature_mean + feature_sd),
            lower_bound_null = pmax(0, null_feature_mean - null_feature_sd),
            upper_bound_null = pmin(1, null_feature_mean + null_feature_sd)
        )
    } else {
        df = df |> mutate(
            lower_bound = feature_mean - feature_sd,
            upper_bound = feature_mean + feature_sd,
            lower_bound_null = null_feature_mean - null_feature_sd,
            upper_bound_null = null_feature_mean + null_feature_sd
        )
    }
    
    return(df)
}




