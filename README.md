# Ethinic-Disparity-in-Melanoma

cramer_20210916.R：heatmap of collinearity among categorical covariates measured by Cramer’s V method.

survival_20210917.R：generated KM curve of overall survival (A) Overall survial and (B) melanoma-specific survival curve by race/ethnicity for patients in SEER from 2004-2014.  and PHCoxPH summary table.

rf_20210916.R and gbm_20210916.R：sampled and trained data using RF and GBM method for 50 repetitions; generated ROC, AUC, ACC, and varaible importance for each repetition.iter_50times.RData;gbm_50iter.RData.

importance_COMB_20210916.R and mean_roc_func_20210916.R: average results of variable importance and ROC curve of 50 repetitions using RF and GBM.

melanoma_logistic_20210920.R：generated interaction plot. Odds ratio of death for other races versus NHWs (A)With/without ulceration ;(B) At I/II stage and III/IV stage

melanoma_compare_20210917.R：Distribution of tumor thickness, ulceration, stage and poverty proportion factors by race/ethnicity. (A), Tumot thickness boxplot; (B), Poverty proportion boxplot; (C), Ulceration barplot; (D), Stage barplot.

meta_stage.R：Forest plots of results of random effects meta-analysis on stage III/IV proportion of Asian and European patients at diagnosis. 
