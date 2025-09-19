plt_roc_curve = function (df_roc, auroc_value = NULL) {
    #' plots roc curve
    
    # visualize x as false-positive-rate
    df_roc = df_roc |>
        mutate(fpr = 1-specificity) |>
        arrange(desc(.threshold)) # fpr 0 -> 1
    
    p = ggplot(df_roc, aes(x=fpr, y=sensitivity)) + 
        geom_step(direction = "vh", color = "black", linewidth = 1) + 
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
        labs(x = "1 - Specificity", y = "Sensitivity") +
        theme_classic(base_size = 16)
    
    if (!is.null(auroc_value)) {
        auroc_value = round(auroc_value, 3)
        p = p + 
            annotate(
                "text",
                label = paste0("AUROC\n", auroc_value),
                x = Inf, y=-Inf, hjust = 1.1, vjust = -0.2,
                size = 8
            )
    }
    
    return(p)
}
