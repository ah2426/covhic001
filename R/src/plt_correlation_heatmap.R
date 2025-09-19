library(ggcorrplot)

plt_correlation_heatmap = function (corr, 
                                    pval=NULL, 
                                    type="full", 
                                    title, 
                                    legend.title="Correlation") {
    # hirarchical clustering
    hc_col = hclust(dist(t(corr)))$order
    corr = corr[hc_col, hc_col]
    pval = pval[hc_col, hc_col]
    
    p = .get_corrplot(corr, pval, type, title=title, legend.title=legend.title)
    return(p)
}


.get_corr = function (df, type="corr.rho") {
    # df: output list object from `d00_innate_corr.R` lst[[ctype]][[group]]
    
    to_rm = setdiff(c("corr.rho", "pval"), type)

    mat = df |> select(-all_of(c(to_rm, "celltype", "group"))) |>
        pivot_wider(names_from="pathway", values_from=!!type) |>
        tibble_to_matrix("atac_pred")

    return(mat)
}


.get_corrplot = function (corr, pval, type="full", title, legend.title="Corr") {
    if (is.null(pval) == FALSE) {
        # pvalue reconfigure **due to ggcorrplot feature
        pval[is.na(pval)] = 2 # also significant
        pval[pval < .05] = 2 # significant
        pval[pval >=.05 & pval <= 1] = 0
        pval[pval == 2] = 1

        pval = t(pval)
    }
    
    p = ggcorrplot(t(corr), p.mat=pval,
                   outline.color="white", type=type, 
                   pch.cex=2, pch.col="black", pch=8, legend.title=legend.title) +  # *
        theme_classic(base_size=12) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+
        labs(x="", y="", title=title)

    return(p)
}
