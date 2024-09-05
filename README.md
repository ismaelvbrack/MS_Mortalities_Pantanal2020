# MS_Mortalities_Pantanal2020
This repository contains data and code scripts for the manuscript "*Spatial modelling and estimation of mammals’ mortalities by Pantanal 2020 megafires*" (Brack et al. 2024, **Journal of Applied Ecology**).

## Repository content

-   `data/`: original data
    - `mammal_info.txt`: information of the mammal species used in the study.
    - `mammal_records.txt`: observed data: recorded detections of mammalian carcasses using the dependent double observer protocol.
    - `plots_covariates.txt`: site-level covariates, measured for each sampling plot (1-ha quadrat) with its respective ID.

-   `figs/`: figures generated from the data analysis
    -   `Fig.Lambda-dBurnAWB.png`: Relationship of the expected number of carcasses with the wildfire severity for sites with and without artificial water body.
    -   `Fig.Psi_Coefs_Forest.png`: Slope coefficient of the effect of non-flooded forests in the suitability of carcasses for the entire community and for the evaluated taxa.
    -   `Fig.Psi-Forest.png`: Relationship of suitability of carcasses with the proportion of non-flooded forests. Community average relationship with 95%CI shadow in green, and mean relationship for each taxa.

-   `ms/`: manuscript files used in the submission
  
-   `outputs/`: output objects from models fitted in JAGS  

-    `R/`: R code scripts
    -   `JAGS_models/`: folder containing JAGS models.
        -   `multisppDA_fatalNmix_covars3_6scales.txt`: selected model: psi(forest) lambda(deltaBurn+tanque) p(.)
        -   `multisppDA_fatalNmix_covars3_6scales.txt`: full model: psi(forest+dist2wiw) lambda(deltaBurn+tanque) p(greenVeg)
    -   `figures_relationship predictions.R`: code to make the figures
    -   `results_coefficient estimates.R`: code to see model results
    -   `run_multispp ZIP Nmix carcass.R`: code to import, arrange, and analyze data using JAGS
