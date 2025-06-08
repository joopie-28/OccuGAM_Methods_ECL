#### Simulation framework for OccuGAMs MS

# original code by Nick Clark
library(mvgam)
library(gratia)
library(marginaleffects)
library(dplyr)
library(ggplot2); theme_set(theme_classic())
library(patchwork)
library(scoringRules)


# Source custom functions
source('Functions/012_threshold_functions.R')


#### 1. Create Simulation Parameter Grid

# Define components
species_df <- data.frame(
  species = c('sp1', 'sp2', 'sp3', 'sp4'),
  p_detect = c(0.8, 0.2, 0.2, 0.8),
  mean_abundance = c(100, 100, 40, 40)
)

# Add the non-species inputs
thresholds <- data.frame(threshold = c('piecewise', 'linear', 'nonmonotonic','monotonic', 'cubic','hyperabundance'))
models <- data.frame(models= c('linear','quadratic', 'cubic', 'GAM'),
                     model_scenario = c(1,2,3,4))
site_reps <- data.frame(
  n_sites = c(50),
  n_replicates = c(10),
  scenario = c('large')
)
slurm_scenario <- data.frame('slurm_scenario' = c(1:24))
simrep <- data.frame('simrep' = c(1:100))

# Cross species, thresholds, and site/replicate pairs
study_design <- merge(merge(species_df, thresholds), site_reps)
study_design <- cbind(study_design, design_scenario = c(1:24))

param_grid <- merge(merge(study_design, models), simrep)

# keep environment clean
rm(simrep, site_reps, slurm_scenario, species_df,
   thresholds, models)

# Save file for downstream analyses
saveRDS(param_grid, "/Users/sassen/OccuGAM_Methods_ECL/Analyses/HPC_Packages/06_HPC_Simulations/data/param_grid.rds")

##### Create simulated datasets for model fitting

# collapse down by model fitted, so that we fit the 3 model variations to the same simulated data (per scenario-rep combo)
data_grid <- param_grid |>
  select(-models, -model_scenario) |># Remove model-related columns
  distinct()

simData <- apply(study_design, 1, function(row) {

  rep_list <- lapply(1:100, function(rep){
    print(rep)
    simdat <- sim_nmix(threshold = row['threshold'],
                       n_sites = row['n_sites'],
                       n_replicates = row['n_replicates'],
                       p_detect = as.numeric(row['p_detect']),
                       mean_abundance = as.numeric(row['mean_abundance']))
    return(simdat)
  })
  return(rep_list)
})


# Save the file
saveRDS(simData, "/Users/sassen/OccuGAM_Methods_ECL/Analyses/HPC_Packages/06_HPC_Simulations/data/SimulationInputData.rds")


