#######################################################################
### Non-Linear responses to environmental covariates in SEA Mammals ###
#######################################################################

# Methods comparison paper

# J.M. Sassen 
# Ecological Cascades Lab, University of Queensland
# 12-07-2024

# Main Analysis Part 3b: Fit Occupancy models with polynomials on HPC

#### Activate Packages


#### Source Functions


#### HPC Version of main analysis Script

# Used to run the Bayesian occupancy GAMs on the HPC
#install.packages(c('gratia', 'mgcv', 'tidyverse', 'jagsUI', 'abind', 'unmarked'))
library('mgcv',lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.2')
library('dplyr', lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.1')
library('gratia', lib.loc = '/home/s4633921/R/x86_64-pc-linux-gnu-library/4.2')
library(jagsUI, lib.loc = "/home/s4633921/R/x86_64-pc-linux-gnu-library/4.1")
library(abind)
library(unmarked)

# List the custome functions we build here

# A little custom function to add dates to output files
date.wrap <- function(string){
  paste0(string, "_", Sys.Date(), "_JMS")
}

# Negating the 'in' operator for useful functionality
`%!in%` <- Negate(`%in%`)

# Model fitting
# Proper version for jagsui which allows parallel computing
fitoccu.HPC <- function(sp.list, covariate,params.to.monitor,slurm.model.type){
  
  sp <- names(sp.list[1])
  
  if(covariate %!in% colnames(sp.list[[sp]]@siteCovs)){
    
    stop("The chosen covariate is not present in your data")
  }
  # Set up some key parameters
  # Main params
  nsite <- nrow(sp.list[[sp]]@siteCovs)
  n_survey <- ncol(sp.list[[sp]]@y)
  cov <- sp.list[[sp]]@siteCovs[[covariate]]
  
  
  ## TEST
  
  poly.degree <- match(slurm.model.type, c("Linear", "Quadratic", "Cubic")) + 1
  DisCov <- as.matrix(cbind(1, cov, cov^2, cov^3))

  
  # RE params
  Landscape <- as.numeric(sp.list[[sp]]@siteCovs[['Landscape']])
  nLandscape <- length(unique(Landscape))
  
  Source <- as.numeric(sp.list[[sp]]@siteCovs[['source']])
  nSource <- length(unique(Source))
  
  # year was not part of the original data so we need to consturct that here
  
  Year <- as.numeric(as.factor(stringr::str_extract(sp.list[[sp]]@siteCovs[['survey_id']], 
                                                    stringr::regex("(\\d+)(?!.*\\d)"))))
  nYear <-length(unique(Year))
  
  # Observation params
  y <- sp.list[[sp]]@y
  
  # needs some special attention due to umf alteration
  det.cov.vec <- sp.list[[sp]]@obsCovs$num_cams_active_at_date
  new.columns <- ncol(y)
  det.cov <-matrix(det.cov.vec, ncol = new.columns, byrow = TRUE)
  
  #Jags does not like NAs, we have already standardised so ok to replace with 0, but let's keep NA first
  det.cov[is.na(det.cov)] <- is.numeric(NA)
  #det.cov[is.na(det.cov)] <- 0
  
  # There is a small issue with NA's in the OP covariate, due to "Xishuangbanna" (China) not being included
  # in the CRISP database. Luckily, it can be easily verified that Xishuangbanna has NO OIL PALM. So we will
  # set this to the minimum value of the standardised OP covariate, which corresponds to 0.
  cov[is.na(cov)] <- min(cov, na.rm=T) # Fixed
  
  # nice message for clarity
  print(paste0('Fitting Bayesian Occupancy Models to ', sp, ' for ', n_distinct(Landscape), ' distinct landscapes'))
 
  if(slurm.model.type != 'GAM'){
    # Bind all data together in a list for feeding into JAGS
    data_list <- list(
      y = y,
      DisCov = DisCov[,1:poly.degree],
      nsite = nsite,
      nsurvey = n_survey,
      nLandscape = nLandscape,
      Landscape = Landscape,
      nYear = nYear,
      Year = Year,
      nSource = nSource,
      Source = Source,
      Effort = det.cov
    )
    
    # Declare the initial values function for use in JAGS
    my_inits <- function(chain){
      gen_list <- function(chain = chain){
        list(
          z = rep(1, nsite),
          a0 = rnorm(1),
          aEffort = rnorm(1),
          b = rnorm(poly.degree,mean =0,sd=0.1),
          sd.p = runif(1,0, 0.25),
          bLandscape=rnorm(nLandscape),
          sigma.bLandscape= runif(1, 0, 1),
          bYear=rnorm(nYear),
          sigma.bYear= runif(1, 0, 1),
          aSource=rnorm(nSource),
          sigma.aSource= runif(1, 0, 1),
          .RNG.name = switch(chain,
                             "1" = "base::Wichmann-Hill",
                             "2" = "base::Wichmann-Hill",
                             "3" = "base::Super-Duper",
                             "4" = "base::Mersenne-Twister",
                             "5" = "base::Wichmann-Hill",
                             "6" = "base::Marsaglia-Multicarry",
                             "7" = "base::Super-Duper",
                             "8" = "base::Mersenne-Twister"
          ),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,
                    "1" = gen_list(chain),
                    "2" = gen_list(chain),
                    "3" = gen_list(chain),
                    "4" = gen_list(chain),
                    "5" = gen_list(chain),
                    "6" = gen_list(chain),
                    "7" = gen_list(chain),
                    "8" = gen_list(chain)
      )
      )
    }
    
    #jagsui needs to run this prior to input 
    my_inits.val <- lapply(1:nc, my_inits)
  }
  
  # Change to GAM
  if(slurm.model.type == 'GAM'){
    # In order to create an occupancy GAM, we leverage this nifty jagam function,
    # which provides us with a JAGS template for smooths, which we will later alter
    # to build a GAM that model imperfect detection 
    tmp_jags <- mgcv::jagam(response ~ s(cov, k = 5, bs = 'tp'),
                            data = data.frame(
                              response = rep(1, nsite),
                              cov = cov
                            ),
                            family = "binomial",
                            file = "tmp.jags")
    
    # Bind all data together in a list for feeding into JAGS
    data_list <- list(
      y = y,
      X = tmp_jags$jags.data$X,
      S1 = tmp_jags$jags.data$S1,
      nsite = nsite,
      nsurvey = n_survey,
      nLandscape = nLandscape,
      Landscape = Landscape,
      nYear = nYear,
      Year = Year,
      nSource = nSource,
      Source = Source,
      Effort = det.cov,
      zero = tmp_jags$jags.data$zero
    )
    
    # Declare the initial values function for use in JAGS
    my_inits <- function(chain){
      gen_list <- function(chain = chain){
        list(
          z = rep(1, nsite),
          a0 = rnorm(1),
          aEffort = rnorm(1),
          b = rnorm(
            length(tmp_jags$jags.ini$b),
            tmp_jags$jags.ini$b,
            0.2
          ),
          bLandscape=rnorm(nLandscape),
          sd.p = runif(1,0, 0.25),
          sigma.bLandscape= runif(1, 0, 1),
          bYear=rnorm(nYear),
          sigma.bYear= runif(1, 0, 1),
          aSource=rnorm(nSource),
          sigma.aSource= runif(1, 0, 1),
          lambda = rgamma(2,1,1),
          .RNG.name = switch(chain,
                             "1" = "base::Wichmann-Hill",
                             "2" = "base::Wichmann-Hill",
                             "3" = "base::Super-Duper",
                             "4" = "base::Mersenne-Twister",
                             "5" = "base::Wichmann-Hill",
                             "6" = "base::Marsaglia-Multicarry",
                             "7" = "base::Super-Duper",
                             "8" = "base::Mersenne-Twister"
          ),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,
                    "1" = gen_list(chain),
                    "2" = gen_list(chain),
                    "3" = gen_list(chain),
                    "4" = gen_list(chain),
                    "5" = gen_list(chain),
                    "6" = gen_list(chain),
                    "7" = gen_list(chain),
                    "8" = gen_list(chain)
      )
      )
    }
    
    #jagsui needs to run this prior to input 
    my_inits.val <- lapply(1:nc, my_inits)
  }
  
  
  # fit the model leveraging HPC paraellel computing
  temp.mod <- jagsUI::jags(data = data_list, inits = my_inits.val,
                           n.chains = nc, n.adapt = na, n.iter=ni,parallel = T, 
                           n.thin=nt, n.burnin = nb, parameters.to.save = params.to.monitor,
                           model.file = slurm.model)
  
  jagsmod_samples <- temp.mod$samples 
  
  # Calculate LOOs
  mcmcMat = as.matrix(jagsmod_samples,chains=TRUE)
  mcmc_loglik = mcmcMat[,grep("^log_lik",colnames(mcmcMat))]
  

  
  # Save it! LogLiks are for post-processing
  model.outputs <- list('Species' = sp,
                        'Data_List' = data_list,
                        'Full_model' = temp.mod,
                        'LogLik' = mcmc_loglik,
                        'poly' = poly.degree)
  return(model.outputs)
  
} 

# Relevant extractions
DrawPosteriorJAGS.HPC <- function(JAGS_model, sp.list, covariate){
  
  # Extract species name from the model output, and the vanilla unmarked model
  species <- JAGS_model$Species
  
  # extract the relevant covariate (jags models do not store these
  # so we refer back to our umf lists)
  cov <- sp.list[[species]]@siteCovs[[covariate]]
  cov[is.na(cov)] <- min(cov, na.rm=T) # need to make sure we re-account for NAs as not carried over from model fitting.
  
  # Create new data for Polynomials
  newdat <- data.frame('Intercept' = rep(1,1000),
                       'DisCov' =  seq(min(cov), max(cov), length.out = 1000)) |>
    mutate('DisCov2' = DisCov^2,
           'DisCov3' = DisCov^3)
  
  
  # Create new data for GAMs
  g1<-gam(response ~ s(cov, k = 5, bs = 'tp'),
          data = data.frame(
            response = rep(1, JAGS_model$Data_List$nsite),
            cov = cov),
          family = "binomial")
  
  # Create new data
  newdat_lp <- predict(g1, newdata = data.frame('cov' = newdat[["DisCov"]]) , type = 'lpmatrix')

  #mcmc extraction process
  mod_mcmc <- do.call('rbind', JAGS_model$Full_model$samples)
  
  # Get Model predictions from the MCMC output
  
  # Will calculate this for GAMs as well but we will not use it, so it is OK
  intercept <- do.call(rbind, replicate(nrow(newdat), t(mod_mcmc[,1]), simplify = FALSE))
  
  # Remember to comment this out properly
  poly.degree <- JAGS_model$poly 
  ifelse(!grepl("GAM", slurm.model, ignore.case = TRUE),
         tmp_est<-plogis(as.matrix(newdat[,1:poly.degree]) %*% t(mod_mcmc[,1:poly.degree])),
         # GAM approach
         tmp_est <- plogis((newdat_lp %*% t(mod_mcmc[,1:5]))))
  
  # Make sure we get the full intervals
  tmp_est_df <- as.data.frame(t(apply(tmp_est, 1, 
                                      quantile, 
                                      probs = c(0.025,0.5,0.975))))
  
  # Add in our smooth covariate
  tmp_est_df$cov <- newdat$DisCov
  
  # Order if we want to do plotting later
  tmp_est_df<-tmp_est_df[order(tmp_est_df$cov),]
  
  # bundle up
  output <- list('Draws' = list('Quantiles'= tmp_est_df),
                 'Species' = species)
  return(output)
  
}

###------ Set up HPC parameters

#### Read external values from SLURM into R
setting = Sys.getenv("SETTING")  # MCMC setting

# MCMC settings, based on assignment above
## Want burn-in to be ~20% of iterations and then thin = (ni - nb) / ideal n.eff (per chain), ideally 30000 in the long one. 
### Assess n.eff via (ni - nb)/nt * nc 
if(setting == "SHORT"){
  ni <- 3000;  nt <- 5; nb <- 60; nc <- 4; na = NULL      #quick test to make sure code works
}

if(setting == "THESIS"){
  ni = 20000;  nt = 10; nb = 5000 ; nc <- 4; na = NULL   #examine parameter values -
}

if(setting == "MIDDLE"){
  ni = 50000;  nt = 20; nb = 10000 ; nc <- 4; na = NULL   #examine parameter values -
}
if(setting == "HALF_LONG"){
  ni = 250000;  nt = 20; nb = 50000 ; nc <- 4; na = NULL  #publication quality run 
}

## specify the number of cores to be uses, should be the same as the model 
options(mc.cores = 4)

#### Read Job Array index value into R
slurm = Sys.getenv("SLURM_ARRAY_TASK_ID")
slurm = as.numeric(slurm) #imports as character var, not numeric

# Need to read in the new list
SpCovDF <- read.csv("data/SpCovDF.csv")

# Select the relevant species-covariate-model combination for this task
sp <- SpCovDF[slurm, 'species']
covariate <- SpCovDF[slurm, 'covariate']
slurm.model.type <-SpCovDF[slurm, 'model']  #'GAM', 'LINEAR', ETC
slurm.model <- paste0('code/Models_Occu',slurm.model.type,'.R')  

# Parameter list to save for models
params.to.monitor <- c('b' , 'a0', 'aEffort', 'bLandscape', 'var.bLandscape','bYear','SSEobs', 'SSEsim', 'p.val', 'c.hat','log_lik')

                       
# read in umf.list
umf.list.occu <- readRDS("data/umflistocc.rds")

# Thin the UMF list based on chosen combination
sp.list <- umf.list.occu[sp]

##----- Fit the Bayesian Occupancy Model in JAGS
# model template is pulled from adjacent file

start = Sys.time()

# Model fitting
mod <- fitoccu.HPC(sp.list, covariate, params.to.monitor,slurm.model.type)

# Extract all the necessary parameters
mod.extract <- DrawPosteriorJAGS.HPC(mod, sp.list, covariate)

end = Sys.time()

print(paste("Finished running occupancy GAMM model and derivative computation for: ", sp," - ",covariate, " at ", Sys.time(),
            ". This model took ", round(difftime(end, start, units = "hours"), 4), " hours to be completed.", sep = ""))

# Create a useful naming scheme
modname <- date.wrap(gsub(" ","_",paste0(sp,"_", covariate, "_", slurm.model.type)))

### Now a small module to save relevant model params in a file rather than exporting the whole model
# this dynamically updates the df using slurm index
threshold_mod_summary <- readRDS('data/occupancy_mod_summary.rds') 

# Just filling in relevant params and diagnsotics
threshold_mod_summary[slurm, 'Species']  <- sp
threshold_mod_summary[slurm, 'Covariate'] <- covariate
threshold_mod_summary[slurm, 'Model_Type'] <- slurm.model.type
threshold_mod_summary[slurm, 'DIC'] <- mod$Full_model$DIC
threshold_mod_summary[slurm, 'Bayes_P'] <- mod$Full_model$mean$p.val
threshold_mod_summary[slurm, 'C_Hat'] <- mod$Full_model$mean$c.hat
threshold_mod_summary[slurm, 'modname'] <- modname

# This should give us all we need in < 5mb!

# We will collect these 1 row dataframes and then just bind together on the local drive
saveRDS(threshold_mod_summary, paste0('results/OCCU_RDS/', covariate ,"/", modname,'_mod_summary.rds'))

### Saving the logliks for post-processing
saveRDS(mod$LogLik, paste0('results/LogLiks/', covariate ,"/", modname,'_LogLiks.rds'))

# Plotting output

pathname.D <- paste0("results/PLOTTING_RDS/",covariate,"/", modname, ".rds")

saveRDS(mod.extract, pathname.D)

# End of OCCU HPC Script


