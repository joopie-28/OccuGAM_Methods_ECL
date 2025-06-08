##########################
### Analyse HPC output ###
##########################

# GAMs for Wildlife
# 10/11/2024
# J.M. Sassen

# Packages
library(data.table)
library(hydroGOF)
library(plotrix)
library(tidyverse)
library(loo)
library(scales)


# A little plotting utility function

auto.legend.pos <- function(x,y,xlim=NULL,ylim=NULL) {
  if (dev.cur() > 1) {
    p <- par('usr')
    if (is.null(xlim)) xlim <- p[1:2]
    if (is.null(ylim)) ylim <- p[3:4]
  } else {
    if (is.null(xlim)) xlim <- range(x, na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(y, na.rm = TRUE)
  }
  countIt <- function(a) {
    tl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y >= ylim[1]*a+ylim[2]*(1-a))
    tr <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y >= ylim[1]*a+ylim[2]*(1-a))
    bl <- sum(x <= xlim[1]*(1-a)+xlim[2]*a & y <= ylim[1]*(1-a)+ylim[2]*a)
    br <- sum(x >= xlim[1]*a+xlim[2]*(1-a) & y <= ylim[1]*(1-a)+ylim[2]*a)
    c(topleft=tl,topright=tr,bottomleft=bl,bottomright=br)
  }
  for (k in seq(0.5,0.1,by=-0.05)) {
    a <- countIt(k)
    if (sum(a==0)>0) break
  }
  names(a)[which(a==0)][1]   # may delete "[1]"
}

# Specify the results folder
result.folder <- "/Users/sassen/Dropbox/Joop MSc QBIO project/Joop OccuGAMs methods paper/OccuGAMs Methods - Results/Wildlife Case Sudies/05_HPC_Comprehensive" # change to dropbox

# Load in covariates
covariates <- c('Avg_FLLI_3km', 'Avg_forest_cover_3km',
                'Avg_human_footprint_3km', 'Avg_OP_percent_3km')

#############################
### P1 - N-Mixture Models ###
#############################

# Initialise collection list
abu.mod.summary.list <- list()

##### Model Fitting Params ####
# Merge model parameter extracts 

# Loop through the folders to extract our model info, its broken up 
# into different 1-row df's so we need to extract and bind.
for(cov in covariates){
  current.folder <- list.files(paste0(result.folder, "/Abundance/results/ABU_RDS/", cov)) 

  for(mod in current.folder){
    print(mod)
    temp <- readRDS(paste0(result.folder, "/Abundance/results/ABU_RDS/", cov, "/", mod)) 
    slurm.index <- which(temp$Species != "")
    abu.mod.summary.list[[slurm.index]] <- temp[slurm.index,]
  }
}

# Final Modelling Outputs (DICs, P-values, C-hats)
full.abu.mod.summary <- rbindlist(abu.mod.summary.list)

##### NRMSE and Correlation #####
full.abu.mod.summary$correlation_GAM <- NA
full.abu.mod.summary$nrmse_GAM <- NA

# Same Looping functionality but different folder
for(cov in covariates){
  current.folder <- list.files(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov)) 
  
  for(mod in current.folder){
    print(mod)
    temp <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", mod)) 
    
    # find the right entry
    slurm.index <- which(full.abu.mod.summary$Species == temp$Species & full.abu.mod.summary$Covariate == cov & unname(sapply(full.abu.mod.summary$Model_Type, grepl, x = mod, ignore.case = TRUE)))
    
    gc()
    
    # Only do this for non-GAMs
    if(full.abu.mod.summary$Model_Type[slurm.index] != 'GAM'){
      
      # Extract the correct GAM
      GAM.mod <- full.abu.mod.summary$modname[which(full.abu.mod.summary$Species == temp$Species & full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Model_Type == 'GAM')]
      GAM.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", GAM.mod, ".rds"))$Draws$Quantiles
      
      # Extract the correct comparison model
      mod.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", full.abu.mod.summary$modname[slurm.index], ".rds"))$Draws$Quantiles
      
      # NRMSE - normalised difference metric
      
      full.abu.mod.summary$nrmse_GAM[slurm.index] <- nrmse(mod.output$`50%`, GAM.output$`50%`, norm = "maxmin") 

    }
  }
}

# add labels
full.abu.mod.summary$cov_name <- ifelse(full.abu.mod.summary$Covariate == "Avg_human_population_3km", 'Human Population',
                                        ifelse(full.abu.mod.summary$Covariate == 'Avg_OP_percent_3km', 'Oil Palm', 
                                               ifelse(full.abu.mod.summary$Covariate == "Avg_forest_cover_3km", 'Forest Cover', 
                                                      ifelse(full.abu.mod.summary$Covariate == "Avg_human_footprint_3km", 'Human Footprint', 
                                                             'Forest Integrity'))))
# Keep environemnt clean
rm(abu.mod.summary.list, GAM.output, mod.output, temp,
   cov, current.folder, GAM.mod, mod, slurm.index)

# Save the results dataframe
saveRDS(full.abu.mod.summary, "Outputs/Main Results/NMixture_Output.rds")
write_csv(full.abu.mod.summary, "Outputs/Main Results/NMixture_Output.csv")

#############################
### P2 - Occupancy Models ###
#############################

# Initialise collection list
occu.mod.summary.list <- list()

##### Model Fitting Params ####
# Merge model parameter extracts 

# Loop through the folders to extract our model info, its broken up 
# into different 1-row df's so we need to extract and bind.
for(cov in covariates){
  current.folder <- list.files(paste0(result.folder, "/Occupancy/results/OCCU_RDS/", cov)) 
  
  for(mod in current.folder){
    print(mod)
    temp <- readRDS(paste0(result.folder, "/Occupancy/results/OCCU_RDS/", cov, "/", mod)) 
    slurm.index <- which(temp$Species != "")
    occu.mod.summary.list [[slurm.index]] <- temp[slurm.index,]
  }
}

# Final Modelling Outputs (DICs, P-values, C-hats)
full.occu.mod.summary <- rbindlist(occu.mod.summary.list )

##### NRMSE and Correlation #####
full.occu.mod.summary$correlation_GAM <- NA
full.occu.mod.summary$nrmse_GAM <- NA

# Same Looping functionality but different folder
for(cov in covariates){
  current.folder <- list.files(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov)) 
  
  for(mod in current.folder){
    print(mod)
    temp <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov, "/", mod)) 
    
    # find the right entry
    slurm.index <- which(full.occu.mod.summary$Species == temp$Species & full.occu.mod.summary$Covariate == cov & unname(sapply(full.occu.mod.summary$Model_Type, grepl, x = mod, ignore.case = TRUE)))
    
    gc()
    
    # Only do this for non-GAMs
    if(full.occu.mod.summary$Model_Type[slurm.index] != 'GAM'){
      
      # Extract the correct GAM
      GAM.mod <- full.occu.mod.summary$modname[which(full.occu.mod.summary$Species == temp$Species & full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Model_Type == 'GAM')]
      GAM.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov, "/", GAM.mod, ".rds"))$Draws$Quantiles
      
      # Extract the correct comparison model
      mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov, "/", full.occu.mod.summary$modname[slurm.index], ".rds"))$Draws$Quantiles
      
      # NRMSE - normalised difference metric
      
      full.occu.mod.summary$nrmse_GAM[slurm.index] <- nrmse(mod.output$`50%`, GAM.output$`50%`, norm = "maxmin") 
      
    }
  }
}

# add labels
full.occu.mod.summary$cov_name <- ifelse(full.occu.mod.summary$Covariate == "Avg_human_population_3km", 'Human Population',
                                        ifelse(full.occu.mod.summary$Covariate == 'Avg_OP_percent_3km', 'Oil Palm', 
                                               ifelse(full.occu.mod.summary$Covariate == "Avg_forest_cover_3km", 'Forest Cover', 
                                                      ifelse(full.occu.mod.summary$Covariate == "Avg_human_footprint_3km", 'Human Footprint', 
                                                             'Forest Integrity'))))

# Keep environment clean
rm(occu.mod.summary.list, GAM.output, mod.output, temp,
   cov, current.folder, GAM.mod, mod, slurm.index)

# Save the results dataframe
saveRDS(full.occu.mod.summary, "Outputs/Main Results/Occupancy_Output.rds")
write_csv(full.occu.mod.summary, "Outputs/Main Results/Occupancy_Output.csv")

# Export statistics for both Occupancy and Rel. Abundance in a neat table format

for (sp in c('Sus scrofa', 'Macaca nemestrina', 'Rusa unicolor', 'Muntiacus')){
  for(cov in covariates ){
    print(sp)
    print(cov)
   
  # Occupancy
   temp <- full.occu.mod.summary |>
      filter(Species == sp & Covariate == cov) |> 
      mutate_if(is.numeric, round, digits =2) |>
      select(Species,
             Covariate,
             Model_Type,
             nrmse_GAM,
             Bayes_P,
             C_Hat)
   
   write_csv(temp,
             paste0('Outputs/Model_metrics/Occupancy/', sp, '_', cov, '.csv'))
   
   # Abundance
   temp <- full.abu.mod.summary |>
     filter(Species == sp & Covariate == cov) |>
     mutate_if(is.numeric, round, digits =2) |>
     select(Species,
            Covariate,
            Model_Type,
            nrmse_GAM,
            Bayes_P,
            C_Hat)
   
   write_csv(temp,
             paste0('Outputs/Model_metrics/Abundance/', sp, '_', cov, '.csv'))
    
  }
}

# Keep env clean
rm(temp, sp, cov)

##########################
#### Visualise Results ###
##########################

# Step 0. Read species standardisation parameters

sp_standardisation_params <- readRDS("/Users/sassen/OccuGAM_Methods_ECL/Inputs/Plotting/standardisation_params.rds") # note, this is for the reduced landscapes in muntjak

############################
########## Main Figures
############################

#### Figure 1.
# A side-by-side comparison of GAMs and Polynomials
# Oil Palm was chosen for the main results as prior work has shown that some species
# may achieve hyperabundance in OP plantations. It also a very topical environmental factor
# in the context of SEA.

PlotAllModelsCurves <- function(sp, cov, conf.interval =T, comp = 'cubic'){
  cov.name <- ifelse(cov == "Avg_human_population_3km", 'Human Population',
                     ifelse(cov == 'Avg_OP_percent_3km', 'Oil Palm Percent', 
                            ifelse(cov == "Avg_forest_cover_3km", 'Forest Cover (%)', 
                                   ifelse(cov == "Avg_human_footprint_3km", 'Human Footprint (%)', 'Forest Integrity'))))
  
  # revert standardisation for plotting
  mean <- sp_standardisation_params[[sp]][sp_standardisation_params[[sp]]$Variable == cov,2]
  sd <- sp_standardisation_params[[sp]][sp_standardisation_params[[sp]]$Variable == cov,3]
  
  ### Occupancy
  y.max <- 1
  
  layout(matrix(c(1,2,3), ncol = 3, byrow = TRUE), widths = c(1,1, 1,1), heights = c(1.5, 1.5, 1.5, 2.25,2.25))
  par(mar = c(6, 7, 4, 1)) 
  
  # Extract all model types
  gam.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'GAM')
  linear.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Linear')
  quad.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Quadratic')
  cubic.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Cubic')
  
  GAM.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[gam.i], "/", full.occu.mod.summary$modname[gam.i], ".rds"))$Draws$Quantiles
  lin.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[linear.i], "/", full.occu.mod.summary$modname[linear.i], ".rds"))$Draws$Quantiles
  quad.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[quad.i], "/", full.occu.mod.summary$modname[quad.i], ".rds"))$Draws$Quantiles
  cubic.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[cubic.i], "/", full.occu.mod.summary$modname[cubic.i], ".rds"))$Draws$Quantiles
  
  # back transform
  GAM.output$cov <- ((GAM.output$cov * sd) + mean)
  lin.output$cov <- ((lin.output$cov * sd) + mean)
  quad.output$cov <- ((quad.output$cov * sd) + mean)
  cubic.output$cov <- ((cubic.output$cov * sd) + mean)
  
  plot(GAM.output$`50%`~GAM.output$cov, type='n', ylim=c(0, y.max), main = 'Occupancy',ylab = 'Occupancy', xlab = cov.name,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2.5,    # Increase axis label size
       cex.axis = 2.5,
       yaxs="i",
       las=1,
       mgp = c(4.85,1.4,0))   # Increase axis tick label size))
  
  if(conf.interval){
    polygon(c(lin.output$cov, rev(lin.output$cov)), c(lin.output$`97.5%`, rev(lin.output$`2.5%`)), col=alpha('green3',0.1), border =NA)
    polygon(c(quad.output$cov, rev(quad.output$cov)), c(quad.output$`97.5%`, rev(quad.output$`2.5%`)), col=alpha('steelblue',0.1), border =NA)
    polygon(c(cubic.output$cov, rev(cubic.output$cov)), c(cubic.output$`97.5%`, rev(cubic.output$`2.5%`)), col=alpha('purple',0.1), border =NA)
    polygon(c(GAM.output$cov, rev(GAM.output$cov)), c(GAM.output$`97.5%`, rev(GAM.output$`2.5%`)), col=alpha('black',0.1), border =NA)
  }
  
  
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(lin.output$`50%`~lin.output$cov, lty =1,col = 'green3', lwd=2)
  lines(quad.output$`50%`~quad.output$cov, lty =1,col = 'steelblue', lwd=2)
  lines(cubic.output$`50%`~cubic.output$cov, lty =1,col = 'purple', lwd=2)
  
  
  pos <- auto.legend.pos(y=GAM.output$`50%`, x= GAM.output$cov)
  legend(pos, legend = c('GAM','Linear', 'Quadratic', 'Cubic'), lty = c(1, 1,1,1), col=c('black','green3', 'steelblue','purple'), cex=1.5)
  
  ### Abundance
  
  # Extract all model types
  gam.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'GAM')
  linear.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'Linear')
  quad.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'Quadratic')
  cubic.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp &full.abu.mod.summary$Model_Type == 'Cubic')
  
  GAM.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[gam.i], "/", full.abu.mod.summary$modname[gam.i], ".rds"))$Draws$Quantiles
  lin.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[linear.i], "/", full.abu.mod.summary$modname[linear.i], ".rds"))$Draws$Quantiles
  quad.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[quad.i], "/", full.abu.mod.summary$modname[quad.i], ".rds"))$Draws$Quantiles
  cubic.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[cubic.i], "/",full.abu.mod.summary$modname[cubic.i], ".rds"))$Draws$Quantiles
  
  y.max <- ifelse(sp == 'Sus scrofa' | sp == 'Macaca nemestrina', 200, max(GAM.output$`97.5%`) * 1.2 )
  
  # back transform
  GAM.output$cov <- ((GAM.output$cov * sd) + mean)
  lin.output$cov <- ((lin.output$cov * sd) + mean)
  quad.output$cov <- ((quad.output$cov * sd) + mean)
  cubic.output$cov <- ((cubic.output$cov * sd) + mean)
  
  plot(GAM.output$`50%`~GAM.output$cov, type='n', ylim=c(0, y.max), main = 'Abundance',ylab = 'Relative Abundace', xlab = cov.name,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2.5,    # Increase axis label size
       cex.axis = 2.5,
       yaxs="i",
       las=1,
       mgp = c(4.85,1.4,0))   # Increase axis tick label size))
  
  if(conf.interval){
    polygon(c(lin.output$cov, rev(lin.output$cov)), c(lin.output$`97.5%`, rev(lin.output$`2.5%`)), col=alpha('green3',0.1), border =NA)
    polygon(c(quad.output$cov, rev(quad.output$cov)), c(quad.output$`97.5%`, rev(quad.output$`2.5%`)), col=alpha('steelblue',0.1), border =NA)
    polygon(c(cubic.output$cov, rev(cubic.output$cov)), c(cubic.output$`97.5%`, rev(cubic.output$`2.5%`)), col=alpha('purple',0.1), border =NA)
    polygon(c(GAM.output$cov, rev(GAM.output$cov)), c(GAM.output$`97.5%`, rev(GAM.output$`2.5%`)), col=alpha('black',0.1), border =NA)
  }
  
  
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(lin.output$`50%`~lin.output$cov, lty =1,col = 'green3', lwd=2)
  lines(quad.output$`50%`~quad.output$cov, lty =1,col = 'steelblue', lwd=2)
  lines(cubic.output$`50%`~cubic.output$cov, lty =1,col = 'purple', lwd=2)
  
  
  pos <- auto.legend.pos(y=GAM.output$`50%`, x= GAM.output$cov)
  legend(pos, legend = c('GAM','Linear', 'Quadratic', 'Cubic'), lty = c(1, 1,1,1), col=c('black','green3', 'steelblue','purple'), cex=1.5)
  
  # Comparison panel
  if(comp == 'Linear'){
    mod.output <- lin.output
    mod.col <- 'green3'
  } else if( comp == 'Quadratic'){
    mod.output <- quad.output
    mod.col <- 'steelblue'
  } else if( comp == 'Cubic'){
    mod.output <- cubic.output
    mod.col <- 'purple'
  }
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0, y.max ), main = paste0(comp, ' versus GAM'),ylab = 'Relative Abundance', xlab = cov.name,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2.5,    # Increase axis label size
       cex.axis = 2.5,
       yaxs="i",
       las=1,
       mgp = c(4.85,1.4,0))   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), 
          c(mod.output$`50%`, rev(GAM.output$`50%`)), col=alpha('grey',0.4), border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = mod.col, lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c(comp, 'GAM', '95%-CI'), lty = c(1, 1,2), col=c(mod.col,'black', 'grey'), cex=1.5)
  
}

for(sp in c('Sus scrofa', 'Macaca nemestrina', 'Rusa unicolor', 'Muntiacus')){
  for(cov in c('Avg_OP_percent_3km')){
    pdf(paste0('Outputs/Abundance_Plots/', sp, '_',cov, '.pdf'), width =12.5, height=5)
    PlotAllModelsCurves(sp, cov, comp='Linear')
    dev.off()
  }
}

#### Figure 2.
# Predictive performance as calculated by LOO (Leave-one-out) cross validation
# for all abundance models

plotLOOEstimates.ABU <- function(sp){
  
  print(paste0('Commencing leave-one-out CV for ', sp))
  
  res <- list()
  for(cov in c('Avg_FLLI_3km', 'Avg_forest_cover_3km',
               'Avg_human_footprint_3km', 'Avg_OP_percent_3km')){
    current.folder <- list.files(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov)) 
    current.folder <- current.folder[grepl(sp,current.folder, ignore.case = TRUE)]
    
    for(mod in current.folder){
      print(mod)
      temp <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", mod)) 
      
      # find the right entry
      slurm.index <- which(full.abu.mod.summary$Species == temp$Species & full.abu.mod.summary$Covariate == cov & unname(sapply(full.abu.mod.summary$Model_Type, grepl, x = mod, ignore.case = TRUE)))
      
      ## Use log likelihoods to calculate LOOIC
      LLs <- readRDS(paste0(result.folder, "/Abundance/results/LogLiks/", cov, "/", full.abu.mod.summary$modname[slurm.index], "_LogLiks.rds"))
      
      # run LOO-CV
      res[[paste0(temp$Species,"_", cov)]][[mod]] <- loo::loo(LLs) 
      
      
      # Memory management
      rm(LLs)
      gc()
    }
  }
  print('Completed leave-one-out CV')
  
  res2<-lapply(res, FUN = loo_compare)
  
  ######
  q<-lapply(res2, FUN = function(pair) {
    
    
    # Extract the corresponding data frame from res2 using pair, and select the first two columns
    t <- as.data.frame.matrix(pair)[, 1:2]  # Use [[ to correctly reference list elements
    t$name <- rownames(t)
    t$ratio <- t$elpd_diff/t$se_diff
    return(t)
  })
  
  final <- rbindlist(q)
  
  library(stringr)
  
  # Input string
  
  
  # Define possible options for each column
  species_options <- c('Sus_scrofa', 'Macaca_nemestrina', 'Rusa_unicolor', 'Muntiacus')
  habitat_options <- c('Avg_FLLI_3km', 'Avg_forest_cover_3km', 'Avg_human_footprint_3km', 'Avg_OP_percent_3km')
  model_options <- c('GAM', 'Linear', 'Quadratic', 'Cubic')
  
  for(i in 1:nrow(final)){
    # Match the components using str_detect
    final$species[i] <- species_options[str_detect(final$name[i], species_options)]  # Match species
    final$cov[i] <- habitat_options[str_detect(final$name[i], habitat_options)]  # Match habitat
    final$model[i] <- model_options[str_detect(final$name[i], model_options)]  # Match model
  }
  
  final$cov_name <- ifelse(final$cov == "Avg_human_population_3km", 'Human Population',
                           ifelse(final$cov == 'Avg_OP_percent_3km', 'Oil Palm', 
                                  ifelse(final$cov == "Avg_forest_cover_3km", 'Forest Cover', 
                                         ifelse(final$cov == "Avg_human_footprint_3km", 'Human Footprint', 
                                                'Forest Integrity'))))
  
  data <- final
  # change NAs to 0
  # Define the mapping for categories to numbers
  category_mapping <- c("Forest Integrity" = 4, "Forest Cover" = 3, "Human Footprint" = 2, "Oil Palm" = 1)
  color_mapping <- c("GAM" = 'black', "Linear" = 'green3', "Quadratic" = 'steelblue', "Cubic" = 'purple')
  
  # Add a new column 'category_number' based on the text values in 'category'
  data$category_number <- category_mapping[data$cov_name]
  data$color <- color_mapping[data$model]
  
  data[is.na(data)] <-0
  
  data <- data |> 
    arrange(desc(category_number), desc(ratio)) |> 
    group_by(category_number, species) |> 
    mutate(jittered_y = category_number + (row_number() - mean(row_number())) * 0.1) |> 
    ungroup()
  
  
  
  plot.data <- data[data$species == sp,]
  plot.data$ratio <- plot.data$ratio * -1
  
  # Set up the base R plot
  par(mar = c(5.1,15.1,4.1,2.1))
  plot(NA, xlim = c(0,10), ylim = c(0.5, 4.5),
       xlab = expression(Delta ~ "elpd"['loo'] ~ "(" * sigma * ")"), ylab = ""
       , xaxt = "n", yaxt = "n", bty = "n",
       cex.lab=2)
  
  # Add the x-axis with custom ticks
  axis(1,
       cex.axis = 2)
  
  # Add the y-axis with tick labels
  axis(2, at = 4:1, labels = c("Forest Integrity", "Forest Cover", "Human Footprint", "Oil Palm"), las = 1,cex.axis = 2)
  
  # Add vertical lines from 0.1 to the points, with jittered y-values
  segments(x0 = 0, y0 = plot.data$jittered_y, 
           x1 = ifelse(plot.data$ratio-0.1 < 0,plot.data$ratio, plot.data$ratio-0.1 ), y1 = plot.data$jittered_y, 
           col = alpha(plot.data$color, ifelse(plot.data$ratio < 2, 1, 0.4)), lwd = 1.5)
  
  # Plot the models (dots) with jittered y-values
  points(plot.data$jittered_y~ (plot.data$ratio), pch = 21, col = 'black', 
         bg = alpha(plot.data$color, ifelse(plot.data$ratio < 2, 1, 0.5)), cex = 1.5)
  
  box(col='black')
  
  abline(v=2, lwd =1, lty=2, col=alpha('black', 0.4))
  
  
}

for(sp in c('Sus_scrofa', 'Macaca_nemestrina', 'Rusa_unicolor', 'Muntiacus')){
  pdf(paste0("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Loo Performance Plots/", sp, "ABU.pdf"), width = 10, height = 7)
  plotLOOEstimates.ABU(sp)
  dev.off()
}

#### Figure 3.
# The NRMSE (i.e. standardised difference) between GAMs and polynomials. Helps the reader build 
# intuition on how the models relate, without having to plot all curves, for all abundance models

full.abu.mod.summary$Model_Type <- factor(full.abu.mod.summary$Model_Type, levels = c('GAM', 'Linear',
                                                                                      'Quadratic', 'Cubic'))
noGAM<-full.abu.mod.summary[full.abu.mod.summary$Model_Type != 'GAM', ]

pdf('Outputs/Main Figures/Fig2.pdf', width = 12, height =7.5)
ggplot(aes(y=nrmse_GAM, x=fct_relevel(cov_name, "Forest Integrity","Forest Cover","Human Footprint", "Oil Palm"), fill = Model_Type), 
       data=noGAM) + facet_wrap(~fct_relevel(Species, 'Macaca nemestrina', 'Sus scrofa','Rusa unicolor', 'Muntiacus'),
                                             scales="free_y") +
  geom_bar(stat = "identity", 
           position = "dodge",color = "black", size = 0.5)  +
  theme_classic() + scale_fill_manual(values = c("Linear" = "green3", 
                                                 "Quadratic" = "steelblue",
                                                 "Cubic" = 'purple')) +
  xlab('Covariate') + labs(fill = 'Model') +
  ylab('Difference from GAM (NRMSE) ') + theme(
    panel.border = element_rect(color = "black", size = 1, fill ='NA'),
    # plot.background = element_rect(color = "black", size = 1),
    strip.background = element_rect(color = "black", size = 1),  # Adds a box around the facet labels
    strip.text = element_text(size = 16),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12))

dev.off()

library(dplyr)
library(ggplot2)
library(forcats)

# Add dummy y-limit values by species
noGAM <- noGAM %>%
  mutate(y_limit = ifelse(Species %in% c('Macaca nemestrina', 'Sus scrofa'), 30, 200))

pdf('Outputs/Main Figures/Fig2.pdf', width = 12, height =7.5)
# Plot with geom_blank to force custom y-limits
ggplot(aes(
  y = nrmse_GAM,
  x = fct_relevel(cov_name, "Forest Integrity", "Forest Cover", "Human Footprint", "Oil Palm"),
  fill = Model_Type
), data = noGAM) +
  facet_wrap(~fct_relevel(Species, 'Macaca nemestrina', 'Sus scrofa', 'Rusa unicolor', 'Muntiacus'),
             scales = "free_y", nrow = 2) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  geom_blank(aes(y = y_limit)) +  # dummy geom to set y-limits
  theme_classic() +
  scale_fill_manual(values = c("Linear" = "green3",
                               "Quadratic" = "steelblue",
                               "Cubic" = 'purple')) +
  xlab('Covariate') +
  labs(fill = 'Model') +
  ylab('Difference from GAM (NRMSE)') +
  theme(
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    strip.background = element_rect(color = "black", size = 1),
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
dev.off()







############################
########## Supplementary Figures
############################

#### Figures S1-S3.
# Copies of Main Figure 1, with remaining coviaraints

for(sp in c('Sus scrofa', 'Macaca nemestrina', 'Rusa unicolor', 'Muntiacus')){
  for(cov in c('Avg_FLLI_3km', 'Avg_forest_cover_3km', 'Avg_human_footprint_3km')){
    pdf(paste0('Outputs/Abundance_Plots/', sp, '_',cov, '.pdf'), width =12.5, height=5)
    PlotAllModelsCurves(sp, cov, comp='Linear')
    dev.off()
  }
}

#### Figure S4
# The NRMSE (i.e. standardised difference) between GAMs and polynomials. Helps the reader build 
# intuition on how the models relate, without having to plot all curves, for all OCCUPANCY models

full.occu.mod.summary$Model_Type <- factor(full.occu.mod.summary$Model_Type, levels = c('GAM', 'Linear',
                                                                                        'Quadratic', 'Cubic'))
noGAM<-full.occu.mod.summary[full.occu.mod.summary$Model_Type != 'GAM' ]

pdf('Outputs/Main Figures/Fig2_occu.pdf', width = 12, height =7.5)
ggplot(aes(y=nrmse_GAM, x=fct_relevel(cov_name, "Forest Integrity","Forest Cover","Human Footprint", "Oil Palm"), fill = Model_Type), 
       data=noGAM) + facet_wrap(~Species) +
  geom_bar(stat = "identity", 
           position = "dodge",color = "black", size = 0.5)  + ylim(0,600) +
  theme_classic() + scale_fill_manual(values = c("Linear" = "green3", 
                                                 "Quadratic" = "steelblue",
                                                 "Cubic" = 'purple')) +
  xlab('Covariate') + labs(fill = 'Model') +
  ylab('Difference from GAM (NRMSE) ') + theme(
    panel.border = element_rect(color = "black", size = 1, fill ='NA'),
    # plot.background = element_rect(color = "black", size = 1),
    strip.background = element_rect(color = "black", size = 1),  # Adds a box around the facet labels
    strip.text = element_text(size = 16),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12))

dev.off()

#### Figure S5
# Predictive performance as calculated by LOO (Leave-one-out) cross validation
# for all OCCUPANCY models

plotLOOEstimates.OCCU <- function(sp){
  
  print(paste0('Commencing leave-one-out CV for ', sp))
  
  res <- list()
  for(cov in c('Avg_FLLI_3km', 'Avg_forest_cover_3km',
               'Avg_human_footprint_3km', 'Avg_OP_percent_3km')){
    current.folder <- list.files(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov)) 
    current.folder <- current.folder[grepl(sp,current.folder, ignore.case = TRUE)]
    
    for(mod in current.folder){
      print(mod)
      temp <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov, "/", mod)) 
      
      # find the right entry
      slurm.index <- which(full.occu.mod.summary$Species == temp$Species & full.occu.mod.summary$Covariate == cov & unname(sapply(full.occu.mod.summary$Model_Type, grepl, x = mod, ignore.case = TRUE)))
      
      ## Use log likelihoods to calculate LOOIC
      LLs <- readRDS(paste0(result.folder, "/Occupancy/results/LogLiks/", cov, "/", full.occu.mod.summary$modname[slurm.index], "_LogLiks.rds"))
      
      # run LOO-CV
      res[[paste0(temp$Species,"_", cov)]][[mod]] <- loo::loo(LLs) 
      
      
      # Memory management
      rm(LLs)
      gc()
    }
  }
  print('Completed leave-one-out CV')
  
  res2<-lapply(res, FUN = loo_compare)
  
  ######
  q<-lapply(res2, FUN = function(pair) {
    
    
    # Extract the corresponding data frame from res2 using pair, and select the first two columns
    t <- as.data.frame.matrix(pair)[, 1:2]  # Use [[ to correctly reference list elements
    t$name <- rownames(t)
    t$ratio <- t$elpd_diff/t$se_diff
    return(t)
  })
  
  final <- rbindlist(q)
  
  library(stringr)
  
  # Input string
  
  
  # Define possible options for each column
  species_options <- c('Sus_scrofa', 'Macaca_nemestrina', 'Rusa_unicolor', 'Muntiacus')
  habitat_options <- c('Avg_FLLI_3km', 'Avg_forest_cover_3km', 'Avg_human_footprint_3km', 'Avg_OP_percent_3km')
  model_options <- c('GAM', 'Linear', 'Quadratic', 'Cubic')
  
  for(i in 1:nrow(final)){
    # Match the components using str_detect
    final$species[i] <- species_options[str_detect(final$name[i], species_options)]  # Match species
    final$cov[i] <- habitat_options[str_detect(final$name[i], habitat_options)]  # Match habitat
    final$model[i] <- model_options[str_detect(final$name[i], model_options)]  # Match model
  }
  
  final$cov_name <- ifelse(final$cov == "Avg_human_population_3km", 'Human Population',
                           ifelse(final$cov == 'Avg_OP_percent_3km', 'Oil Palm', 
                                  ifelse(final$cov == "Avg_forest_cover_3km", 'Forest Cover', 
                                         ifelse(final$cov == "Avg_human_footprint_3km", 'Human Footprint', 
                                                'Forest Integrity'))))
  
  data <- final
  # change NAs to 0
  # Define the mapping for categories to numbers
  category_mapping <- c("Forest Integrity" = 4, "Forest Cover" = 3, "Human Footprint" = 2, "Oil Palm" = 1)
  color_mapping <- c("GAM" = 'black', "Linear" = 'green3', "Quadratic" = 'steelblue', "Cubic" = 'purple')
  
  # Add a new column 'category_number' based on the text values in 'category'
  data$category_number <- category_mapping[data$cov_name]
  data$color <- color_mapping[data$model]
  
  data[is.na(data)] <-0
  
  data <- data |> 
    arrange(desc(category_number), desc(ratio)) |> 
    group_by(category_number, species) |> 
    mutate(jittered_y = category_number + (row_number() - mean(row_number())) * 0.1) |> 
    ungroup()
  
  
  
  plot.data <- data[data$species == sp,]
  plot.data$ratio <- plot.data$ratio * -1
  
  # Set up the base R plot
  par(mar = c(5.1,15.1,4.1,2.1))
  plot(NA, xlim = c(0,10), ylim = c(0.5, 4.5),
       xlab = expression(Delta ~ "elpd"['loo'] ~ "(" * sigma * ")"), ylab = ""
       , xaxt = "n", yaxt = "n", bty = "n",
       cex.lab=2)
  
  # Add the x-axis with custom ticks
  axis(1,
       cex.axis = 2)
  
  # Add the y-axis with tick labels
  axis(2, at = 4:1, labels = c("Forest Integrity", "Forest Cover", "Human Footprint", "Oil Palm"), las = 1,cex.axis = 2)
  
  # Add vertical lines from 0.1 to the points, with jittered y-values
  segments(x0 = 0, y0 = plot.data$jittered_y, 
           x1 = ifelse(plot.data$ratio-0.1 < 0,plot.data$ratio, plot.data$ratio-0.1 ), y1 = plot.data$jittered_y, 
           col = alpha(plot.data$color, ifelse(plot.data$ratio < 2, 1, 0.4)), lwd = 1.5)
  
  # Plot the models (dots) with jittered y-values
  points(plot.data$jittered_y~ (plot.data$ratio), pch = 21, col = 'black', 
         bg = alpha(plot.data$color, ifelse(plot.data$ratio < 2, 1, 0.5)), cex = 1.5)
  
  box(col='black')
  
  abline(v=2, lwd =1, lty=2, col=alpha('black', 0.4))
  
  
}

for(sp in c('Sus_scrofa', 'Macaca_nemestrina', 'Rusa_unicolor', 'Muntiacus')){
  pdf(paste0("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Loo Performance Plots/", sp, "OCCU.pdf"), width = 10, height = 7)
  plotLOOEstimates.OCCU(sp)
  dev.off()
}

#### Figure S6
# Posterior Predictive Checks for Occupancy models

# A - Bayesian P-values
pdf('Outputs/Main Figures/Fig3A_occu.pdf', width = 12, height =7.5)
ggplot(aes(y=Bayes_P, x=fct_relevel(Model_Type, "GAM", "Linear", "Quadratic", "Cubic"), fill = cov_name), 
       data=full.occu.mod.summary) + facet_wrap(~Species) +
  geom_bar(stat = "identity", 
           position = "dodge",color = "black", size = 0.5) +
  ylim(0,1) + ylab('Bayesian P-Value')+ labs(fill = 'Covariate') + xlab('Model Type') +
  theme_classic() + scale_fill_manual(values = c("Forest Integrity" = "darkgreen", 
                                                 "Forest Cover" = "lightgreen", 
                                                 "Human Footprint" = "brown",
                                                 "Oil Palm" = 'yellow4'))+
  geom_hline(yintercept = 0.25, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "black", linetype = "solid") +
  geom_hline(yintercept = 0.75, color = "red", linetype = "dashed")+ theme(
    panel.border = element_rect(color = "black", size = 1, fill ='NA'),
    # plot.background = element_rect(color = "black", size = 1),
    strip.background = element_rect(color = "black", size = 1),  # Adds a box around the facet labels
    strip.text = element_text(size = 16),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12))  # Adds a box around the facet labels

dev.off()

# B - C-hat
pdf('Outputs/Main Figures/Fig3B_occu.pdf', width = 12, height =7.5)
ggplot(aes(y=C_Hat, x=fct_relevel(Model_Type, "GAM", "Linear", "Quadratic", "Cubic"), fill = cov_name), 
       data=full.occu.mod.summary) + facet_wrap(~Species) +
  geom_bar(stat = "identity", 
           position = "dodge",color = "black", size = 0.5) +ylim(0,1.5) + ylab('Overdispersion (C-hat)')+ labs(fill = 'Covariate') + xlab('Model Type') +
  theme_classic() + scale_fill_manual(values = c("Forest Integrity" = "darkgreen", 
                                                 "Forest Cover" = "lightgreen", 
                                                 "Human Footprint" = "brown",
                                                 "Oil Palm" = 'yellow4'))+
  geom_hline(yintercept = 1, color = "black", linetype = "solid") +
  geom_hline(yintercept = 1.1, color = "red", linetype = "dashed") + theme(
    panel.border = element_rect(color = "black", size = 1, fill ='NA'),
    # plot.background = element_rect(color = "black", size = 1),
    strip.background = element_rect(color = "black", size = 1),  # Adds a box around the facet labels
    strip.text = element_text(size = 16),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12))  # Adds a box around the facet labels

dev.off()

##### Table of PPC for occu (3A)

full.occu.mod.summary |>
  mutate(Model = paste0(Species, "~", Covariate ),
         Type = 'Occupancy',
         Bayes_P = round(Bayes_P, 3),
         C_Hat = round(C_Hat, 3)) |>
  select(Model,
         Type,
         Formulation = Model_Type,
         `Bayesian P-value` = Bayes_P,
         `C-Hat` = C_Hat,
         NRMSE = nrmse_GAM) |>
  write.csv("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Main Results/OccuOutputTable.csv")

#### Figure S7
# Posterior Predictive Checks for Abundance models

# A - Bayesian P-values
pdf('Outputs/Main Figures/Fig3A.pdf', width = 12, height =7.5)
ggplot(aes(y=Bayes_P, x=fct_relevel(Model_Type, "GAM", "Linear", "Quadratic", "Cubic"), fill = cov_name), 
       data=full.abu.mod.summary) + facet_wrap(~Species) +
  geom_bar(stat = "identity", 
           position = "dodge",color = "black", size = 0.5) +
  ylim(0,1) + ylab('Bayesian P-Value')+ labs(fill = 'Covariate') + xlab('Model Type') +
  theme_classic() + scale_fill_manual(values = c("Forest Integrity" = "darkgreen", 
                                                 "Forest Cover" = "lightgreen", 
                                                 "Human Footprint" = "brown",
                                                 "Oil Palm" = 'yellow4'))+
  geom_hline(yintercept = 0.25, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 0.5, color = "black", linetype = "solid") +
  geom_hline(yintercept = 0.75, color = "red", linetype = "dashed")+ theme(
    panel.border = element_rect(color = "black", size = 1, fill ='NA'),
    # plot.background = element_rect(color = "black", size = 1),
    strip.background = element_rect(color = "black", size = 1),  # Adds a box around the facet labels
    strip.text = element_text(size = 16),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12))  # Adds a box around the facet labels

dev.off()

##### Table of PPC for Abundance (3A)

full.abu.mod.summary |>
  mutate(Model = paste0(Species, "~", Covariate ),
         Type = 'N-Mixture',
         Bayes_P = round(Bayes_P, 3),
         C_Hat = round(C_Hat, 3)) |>
  select(Model,
         Type,
         Formulation = Model_Type,
         `Bayesian P-value` = Bayes_P,
         `C-Hat` = C_Hat,
         NRMSE = nrmse_GAM) |>
  write.csv("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Main Results/AbundanceOutputTable.csv")


# B - C-hat
pdf('Outputs/Main Figures/Fig3B.pdf', width = 12, height =7.5)
ggplot(aes(y=C_Hat, x=fct_relevel(Model_Type, "GAM", "Linear", "Quadratic", "Cubic"), fill = cov_name), 
       data=full.abu.mod.summary) + facet_wrap(~Species) +
  geom_bar(stat = "identity", 
           position = "dodge",color = "black", size = 0.5) +ylim(0,1.5) + ylab('Overdispersion (C-hat)')+ labs(fill = 'Covariate') + xlab('Model Type') +
  theme_classic() + scale_fill_manual(values = c("Forest Integrity" = "darkgreen", 
                                                 "Forest Cover" = "lightgreen", 
                                                 "Human Footprint" = "brown",
                                                 "Oil Palm" = 'yellow4'))+
  geom_hline(yintercept = 1, color = "black", linetype = "solid") +
  geom_hline(yintercept = 1.1, color = "red", linetype = "dashed") + theme(
    panel.border = element_rect(color = "black", size = 1, fill ='NA'),
    # plot.background = element_rect(color = "black", size = 1),
    strip.background = element_rect(color = "black", size = 1),  # Adds a box around the facet labels
    strip.text = element_text(size = 16),
    axis.title = element_text(size=14),
    axis.text = element_text(size=12))  # Adds a box around the facet labels

dev.off()


#### Legend for LOO Plots
pdf(paste0("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Loo Performance Plots/legend1.pdf"))

# Example setup for the legend
legend_colors <- c('black','green3', 'steelblue', 'purple')    # Replace with your actual colors
legend_labels <- c("GAM", "Linear", "Quadratic",'Cubic') 
legend_points <- c(21, 21, 21, 21)  # Consistent point style

# Create the 2x2 legend
plot.new()  # Start a new blank plot
legend("center", legend = legend_labels, pch = legend_points, 
       col = "black", pt.bg = legend_colors, pt.cex = 2, 
       ncol = 4, bty = "n", title = "")  # Arrange in 2 columns

dev.off()

pdf(paste0("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Loo Performance Plots/legend2.pdf"))
# Create the 2x2 legend

legend_colors <- c('purple', alpha('purple',0.5))    # Replace with your actual colors
legend_labels <- c(expression("< 2"~ sigma),expression("> 2"~ sigma) )
legend_points <- c(21, 21, 21, 21)  # Consistent point style
plot.new()  # Start a new blank plot
legend("center", legend = legend_labels, pch = legend_points, 
       col = "black", pt.bg = legend_colors, pt.cex = 2, 
       ncol = 2, bty = "n", title = "")  # Arrange in 2 columns

dev.off()




