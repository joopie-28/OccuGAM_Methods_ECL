##
##### Run UBMS Models in HPC environment
##

# J.M. Sassen 15-09-2024
unloadNamespace('Rcpp')
library('Rcpp',lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.2')
library('ubms',lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.2')
library('RSpectra',lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.2')

#### Source Functions - this is the HPC environment so we explicitly declare functions in this script

# A little custom function to add dates to output files
date.wrap <- function(string){
  paste0(string, "_", Sys.Date(), "_JMS")
}

# Negating the 'in' operator for useful functionality
`%!in%` <- Negate(`%in%`)

# Set up the data parameters
#### Read Job Array index value into R
slurm = Sys.getenv("SLURM_ARRAY_TASK_ID")
slurm = as.numeric(slurm) #imports as character var, not numeric

# Need to read in the new list
SpCovDF <- read.csv("data/SpCovDF.csv")

# Select the relevant species-covariate pair for this task
sp <- SpCovDF[slurm, 'species']
cov <- SpCovDF[slurm, 'covariate']

# read in umf.list
umf.list.occu <- readRDS("data/umflistocc.rds")

# Fit and plot the polynomial model versions using unmarked with increasing polynomial terms of the covariate
occu.models <- list()
formulae <- list(
  "Linear" = as.formula(paste("~ num_cams_active_at_date ~", cov,"+ (1|Landscape) + (1|Year)")),
  "Quadratic" = as.formula(paste("~ num_cams_active_at_date ~ poly(", cov, ", 2, raw = TRUE) + (1|Landscape) + (1|Year)")),
  "Cubic" = as.formula(paste("~ num_cams_active_at_date ~ poly(", cov, ", 3, raw = TRUE) + (1|Landscape) + (1|Year)")))

### Fit occupancy models
# Set up the year covariate for the data

umf.list.occu[[sp]]@siteCovs$Year <- as.numeric(as.factor(stringr::str_extract(umf.list.occu[[sp]]@siteCovs[['survey_id']], 
                                                                               stringr::regex("(\\d+)(?!.*\\d)"))))
# Occupancy models
for (name in names(formulae)) {
  print(paste0('Fitting frequentist models, currently on the ', name, ' model'))
  
  # Save models in the list and export later
  occu.models[[name]] <- ubms::stan_occu(formulae[[name]], data = umf.list.occu[[sp]])
}

### Save the models in the HPC environment

# Create a useful naming scheme
modname <- date.wrap(gsub(" ","_",paste0(sp,"_Poly_", covariate)))

pathname <- paste0('results/Occupancy/', modname, ".rds")

# Save it on the HPC environment
saveRDS(occu.models, pathname)











