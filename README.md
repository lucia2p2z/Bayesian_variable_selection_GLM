## To run simulations:

To reproduce the results on syntehtic data of the manuscript: 
1. To run: launch the code on files `simulations_with_baselines/config_**_baselines_cd.R` with the desidered hyperparameters.
2. To plot: launch the code `simulations_with_baselines/plot_**_cd.R` with the corresponding path.

with `**` being `binomial`, `gaussian`, `poisson`.


## To run experiments on real data:

To reproduce the results on "Heart Disease" dataset (Logistic regression):
1. To estimate parameters on the whole dataset: run `experiments_heart_data/main_baselines.R`
2. To perform estimation and variable selection with CV: run `experiments_heart_data/main_cv_baselines.R`
3. To evaluate and plot the results: run `experiments_heart_data/compare_results_seed_with_baselines.R`

The same applies for the datasets "Pollution" (Linear regression) and "Crabs" (Poisson regression) dataset. 
