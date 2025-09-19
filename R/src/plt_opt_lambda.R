plt_opt_lambda = function (est, opt_index) {
    df_fit = tibble(sparsity = est$sparsity, lambda = est$lambda)
    df_opt = df_fit[opt_index, ]
    
    p = ggplot(df_fit, aes(x=lambda, y=sparsity)) +
        geom_point(size=3, shape=21, fill="darkgrey", color="white") +
        geom_smooth(method="loess") +
        geom_point(data=df_opt, size=6, color="red", alpha=.8) +
        theme_bw(base_size=12) +
        labs(title = paste0("opt.lambda = ", round(df_opt$lambda,3)))
    
    return(p)
}
