#######################################################################
### Non-Linear responses to environmental covariates in SEA Mammals ###
#######################################################################

# Methods comparison paper

# J.M. Sassen 
# Ecological Cascades Lab, University of Queensland
# 12-07-2024

# Main Analysis Part 2b: Fit Abundance Models on HPC (Unmarked Pcount with ^1 and ^2, and Bayesian Abundance GAMs)

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
fitoccuGAM.HPC <- function(sp.list, covariate){
  
  sp <- names(sp.list[1])
  
  if(covariate %!in% colnames(sp.list[[sp]]@siteCovs)){
    
    stop("The chosen covariate is not present in your data")
  }
  # Set up some key parameters
  # Main params
  nsite <- nrow(sp.list[[sp]]@siteCovs)
  n_survey <- ncol(sp.list[[sp]]@y)
  cov <- sp.list[[sp]]@siteCovs[[covariate]]
  
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
  print(paste0('Fitting the Bayesian Occupancy GAM to ', sp, ' for ', n_distinct(Landscape), ' distinct landscapes'))
  
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
  
  # Load in custome function for EDF
  Calculate.EDF.JAGSUI <- function(temp.mod, tmp_jags){
    
    
    pregam <- tmp_jags$pregam
    test.2 <- list()
    test.2$b <- abind(lapply(temp.mod$samples, function(df) t(df[, 1:5])), along = 3)
    test.2$rho <- abind(lapply(temp.mod$samples, function(df) t(df[, 6:7])), along = 3)
    
    rownames(test.2) <- NULL
    pregam$Vp <- cov(t(test.2$b[,,1]))
    pregam$coefficients <- rowMeans(test.2$b[,,1])
    n.chain <- dim(test.2$b)[3]
    
    if (n.chain>1) { 
      for (i in 2:n.chain) {
        pregam$Vp <- pregam$Vp +  cov(t(test.2$b[,,i]))
        pregam$coefficients <-  pregam$coefficients + rowMeans(test.2$b[,,i])
      }
      pregam$Vp <- pregam$Vp/n.chain
      pregam$coefficients <-  pregam$coefficients/n.chain
    }
    
    XWX <- t(pregam$X) %*% (pregam$X) 
    
    ## tr((X'WX + S)^{-1}X'WX
    rho <- rowMeans(test.2$rho); lambda <- exp(rho)
    XWXS <- XWX
    for (i in 1:length(lambda)) {
      ind <- pregam$off[i]:(pregam$off[i]+ncol(pregam$S[[i]])-1)
      XWXS[ind,ind] <-  XWXS[ind,ind] + pregam$S[[i]] * lambda[i]
    } 
    pregam$edf <-  diag(solve(XWXS,XWX))
    EDF <-  sum(pregam$edf[-1])
    
    return(EDF)
  }
  
  #jagsui needs to run this prior to input 
  my_inits.val <- lapply(1:nc, my_inits)
  
  # fit the model leveraging HPC paraellel computing
  temp.mod <- jagsUI::jags(data = data_list, inits = my_inits.val,
                           n.chains = nc, n.adapt = na, n.iter=ni,parallel = T, 
                           n.thin=nt, n.burnin = nb, parameters.to.save = c('b', 'rho', 'a0', 'aEffort', 'bLandscape', 'var.bLandscape','bYear',
                                                                            'SSEobs', 'SSEsim', 'p.val', 'log_lik'),
                           model.file = "code/01HPC_Models_OccuGAM_vECL.R")
  
  jagsmod_samples <- temp.mod$samples 
  
  # Calculate LOOs
  mcmcMat = as.matrix(jagsmod_samples,chains=TRUE)
  mcmc_loglik = mcmcMat[,grep("^log_lik",colnames(mcmcMat))]

  WAIC<-loo::waic(mcmc_loglik)
  LOO<-loo::loo(mcmc_loglik)
  
  # Calculate EDF using our custom function, which is essentially equivalent to type 0
  # EDF from S. Woods 'Sim2jam' function
  
  edf <- Calculate.EDF.JAGSUI(temp.mod, tmp_jags)
  
  # Create a useful naming scheme
  modname <- date.wrap(gsub(" ","_",paste0(sp,"_", covariate)))
  #pathname <- paste0("Outputs/Occupancy Models/",covariate,"/", modname, ".rds")
  #savreds when hPC file scheme done
  
  # Save it!
  model.outputs <- list('OccuGAM_Model' = jagsmod_samples,
                        'EDF' = edf,
                        'Species' = sp,
                        'Data_List' = data_list,
                        'Full_model' = temp.mod,
                        'WAIC' = WAIC,
                        'LOO' =LOO)
  
  
  
  return(model.outputs)
  
} 

# Relevant extractions
DrawPosteriorJAGS.HPC <- function(JAGS_model, sp.list, covariate){
  
  # Extract EDF and species name from the model output, and the vanilla unmarked model
  EDF <- JAGS_model$EDF
  species <- JAGS_model$Species
  
  
  # extract the relevant covariate (jags models do not store these
  # so we refer back to our umf lists)
  cov <- sp.list[[species]]@siteCovs[[covariate]]
  cov[is.na(cov)] <- min(cov, na.rm=T) # need to make sure we re-account for NAs as not carried over from model fitting.
  
  #mcmc extraction process
  mod_mcmc <- do.call('rbind', JAGS_model$Full_model$samples)
  
  # which landscape has lowest SD?
  # ranefs <- which(grepl("bLandscape",names(JAGS_model$Full_model$summary[,2])))
  # ranefs <- ranefs[-length(ranefs)]
  # n.min.ranef <- which.min(JAGS_model$Full_model$summary[ranefs,'sd'])
  # n.min.df <- ranefs[n.min.ranef]
  
  # Get Model predictions from the MCMC output
  
  tmp_est <- plogis((JAGS_model$Data_List$X %*% t(mod_mcmc[,1:5]))) 
  
  
  # Make sure we get the full intervals
  tmp_est_df <- as.data.frame(t(apply(tmp_est, 1, 
                                      quantile, 
                                      probs = c(0.025,0.5,0.975))))
  
  # Add in our smooth covariate
  tmp_est_df$cov <- cov
  
  # Order if we want to do plotting later
  tmp_est_df<-tmp_est_df[order(tmp_est_df$cov),]
  
  
  
  # Plot some posterior curves, WE ARE NOT DRWING ANYMORE JUST DO ALL!!!!
  randoms <- seq(1, ncol(tmp_est))
  
  # An empty dataframe for our draws
  pos.draws <- data.frame(matrix(NA, nrow = nrow(tmp_est_df),
                                 ncol = length(randoms)))
  
  # Populate the dataframe with n draws from the model
  for(i in 1:length(randoms)){
    
    tmp.curve <- plogis(JAGS_model$Data_List$X %*% mod_mcmc[randoms[i],1:5])
    pos.draws[,i] <- tmp.curve
  }
  
  # add the covariate back in and order dataframe
  pos.draws$cov <- cov
  pos.draws <- pos.draws[order(pos.draws$cov),]
  
  # Skip derivative calculation if linear
  if(EDF <=2){
    return(list('Draws' = list('Quantiles'=tmp_est_df),
                'EDF' = EDF,
                'Species' = species))
  }
  
  # Now save those selected draws, we will use them to calculate derivatives and store them
  
  # instnatiate new dataframe for first and second derivatives
  
  d2.df <- data.frame(matrix(NA, nrow = nrow(pos.draws),
                             ncol = length(randoms)))
  
  for (col in 1:(ncol(d2.df))){
    
    gam.temp <- mgcv::gam(pos.draws[,col]~s(cov, k=-1,bs = 'tp'), data = pos.draws, method ='REML')
    
    # Isolate our preditors for derivative functions (slightly touchy functions)
    newdat =  data.frame('cov'= pos.draws$cov)
    
    # Calculate 1st order derivative
    d2.df[,col] <- derivatives(gam.temp, order = 2, data = newdat,
                               eps=0.00001)[[4]]
    
    
  }
  
  # Add covariate in
  d2.df$cov <- newdat$cov
  
  # Module for extracting Bayesian CIs for each derivative
  deriv.2 <- d2.df[-ncol(d2.df)]
  
  # Line up values and calculate quantiles
  deriv.2.quantiles <- as.data.frame(t(apply(deriv.2, 1, 
                                             quantile, 
                                             probs = c(0.025,0.5,0.975))))
  
  # add the covariate back in
  deriv.2.quantiles$cov <- d2.df$cov
  
  # Run the final GAM for the meidan derivative
  gam.mean <- mgcv::gam(tmp_est_df[,2]~ s(cov, k=-1,bs = 'tp'), 
                        data = tmp_est_df, method ='REML')
  
  # Isolate our preditcors for derivative functions (slightly toughy functions)
  newdat =  data.frame('cov'= tmp_est_df$cov)
  
  # Calculate 2nd order derivative
  d2.mean <- derivatives(gam.mean, order = 2, data = newdat,
                         eps=0.00001)[c(4,9)]
  # Calculate 1st order derivative
  d1.mean <- derivatives(gam.mean, order = 1, data = newdat,
                         eps=0.00001)[c(4,9)]
  
  # bundle up
  output <- list('Draws' = list('Quantiles'=tmp_est_df),
                 'Derivatives' = list('First' = list('Median' = d1.mean),
                                      'Second' =list('Median' = d2.mean, 
                                                     'Quantiles' = deriv.2.quantiles)),
                 'EDF' = EDF,
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

# Select the relevant species-covariate pair for this task
sp <- SpCovDF[slurm, 'species']
covariate <- SpCovDF[slurm, 'covariate']

# read in umf.list
umf.list.occu <- readRDS("data/umflistocc.rds")

# Thin the UMF list based on chosen combination
sp.list <- umf.list.occu[sp]

##----- Fit the Bayesian Occupancy Model in JAGS
# model template is pulled from adjacent file

start = Sys.time()

# Model fitting
mod <- fitoccuGAM.HPC(sp.list, covariate)

# Extract all the necessary parameters
mod.extract <- DrawPosteriorJAGS.HPC(mod, sp.list, covariate)

end = Sys.time()

print(paste("Finished running occupancy GAMM model and derivative computation for: ", sp," - ",covariate, " at ", Sys.time(),
            ". This model took ", round(difftime(end, start, units = "hours"), 4), " hours to be completed.", sep = ""))

# Create a useful naming scheme
modname <- date.wrap(gsub(" ","_",paste0(sp,"_", covariate)))

### Now a small module to save relevant model params in a file rather than exporting the whole model
# this dynamically updates the df using slurm index
threshold_mod_summary <- readRDS('data/occupancy_mod_summary.rds') # this might cause problems in. parallel environemnt


# Just filling in relevant params and diagnsotics
threshold_mod_summary[slurm, 'Species']  <- mod$Species
threshold_mod_summary[slurm, 'Covariate'] <- SpCovDF$covariate[slurm]
threshold_mod_summary[slurm, 'EDF'] <- mod$EDF
threshold_mod_summary[slurm,c('Rhat_B1',
                              'Rhat_B2',
                              'Rhat_B3',
                              'Rhat_B4',
                              'Rhat_B5')] <- mod$Full_model$Rhat$b
threshold_mod_summary[slurm,c('Rhat_A0',
                              'Rhat_AEFFORT')] <- c(mod$Full_model$Rhat$a0, mod$Full_model$Rhat$aEffort)

threshold_mod_summary[slurm,c('Effective_Sample_B1',
                              'Effective_Sample_B2',
                              'Effective_Sample_B3',
                              'Effective_Sample_B4',
                              'Effective_Sample_B5')] <- mod$Full_model$n.eff$b

threshold_mod_summary[slurm,c('Effective_Sample_AO',
                              'Effective_Sample_AEFFORT')] <- c(mod$Full_model$n.eff$a0, mod$Full_model$n.eff$aEffort)

threshold_mod_summary[slurm,c('f_B1',
                              'f_B2',
                              'f_B3',
                              'f_B4',
                              'f_B5')] <- mod$Full_model$f$b

threshold_mod_summary[slurm,c('f_AO',
                              'f_AEFFORT')] <- c(mod$Full_model$f$a0, mod$Full_model$f$aEffort)

threshold_mod_summary[slurm, 'DIC'] <- mod$Full_model$DIC
threshold_mod_summary[slurm, 'Bayes_P'] <- mod$Full_model$mean$p.val

#This should give us all we need in < 5mb!

# We will collect these 1 row dataframes and then just bind together on the local drive
saveRDS(threshold_mod_summary, paste0('results/OCCUGAM_RDS/', covariate ,"/", modname,'_mod_summary.rds'))

### Saving data

# Plotting output

pathname.D <- paste0("results/PLOTTING_RDS/",covariate,"/", modname, ".rds")

saveRDS(mod.extract, pathname.D)

# save full model
# We are also going to save the models, we only run 10 anyways
pathname.E <- paste0("results/FULL_MODELS/", modname, ".rds")
saveRDS(mod, pathname.E)

# End of OCCU HPC Script


