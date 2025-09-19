
# COVHIC001 machine learning analysis

## Procedure

+ data preparation:
    + original data were merged & normalized in `R/0_prep_data/R`
+ model training:
    + elastic-net classifier was trained by running `R/main_elastic_net.R`
    + conditional independence network was estimated by running `R/main_network.R`
+ visualization:
    + visualization of model outputs & related analyses were done using 
      `R/visual_elastic_net.R` and `R/visual_network_ggraph.R`
    + output figures are under `outs/*.pdf`
