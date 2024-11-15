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


# Functions

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
result.folder <- "/Users/sassen/Desktop/05_HPC_Comprehensive" # change to dropbox

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
    
  #  abu.mod.summary.list[[slurm.index]] <- temp[slurm.index,]
    
    # Only do this for non-GAMs
    if(full.abu.mod.summary$Model_Type[slurm.index] != 'GAM'){
      
      # Extract the correct GAM
      GAM.mod <- full.abu.mod.summary$modname[which(full.abu.mod.summary$Species == temp$Species & full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Model_Type == 'GAM')]
      GAM.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", GAM.mod, ".rds"))$Draws$Quantiles
      
      # Extract the correct comparison model
      mod.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", full.abu.mod.summary$modname[slurm.index], ".rds"))$Draws$Quantiles

      # Calculate correlation
  
      full.abu.mod.summary$correlation_GAM[slurm.index] <- cor(mod.output$`50%`, GAM.output$`50%`)
      
      # NRMSE - normalised difference metric
      
      full.abu.mod.summary$nrmse_GAM[slurm.index] <- nrmse(mod.output$`50%`, GAM.output$`50%`, norm = "maxmin") 

    }
  }
}

# Keep environemnt clean
rm(abu.mod.summary.list)


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
    
    
    # Only do this for non-GAMs
    if(full.occu.mod.summary$Model_Type[slurm.index] != 'GAM'){
      
      # Extract the correct GAM
      GAM.mod <- full.occu.mod.summary$modname[which(full.occu.mod.summary$Species == temp$Species & full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Model_Type == 'GAM')]
      GAM.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov, "/", GAM.mod, ".rds"))$Draws$Quantiles
      
      # Extract the correct comparison model
      mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", cov, "/", full.occu.mod.summary$modname[slurm.index], ".rds"))$Draws$Quantiles
      
      # Calculate correlation
      
      full.occu.mod.summary$correlation_GAM[slurm.index] <- cor(mod.output$`50%`, GAM.output$`50%`)
      
      # NRMSE - normalised difference metric
      
      full.occu.mod.summary$nrmse_GAM[slurm.index] <- nrmse(mod.output$`50%`, GAM.output$`50%`, norm = "maxmin") 
      
    }
  }
}

rm(occu.mod.summary.list)

#### Plotting #####

PlotAllModelsOccu <- function(sp, cov){
  
  gam.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'GAM')
  
  GAM.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[gam.i], "/", full.occu.mod.summary$modname[gam.i], ".rds"))$Draws$Quantiles
  
  layout(matrix(c(1,1,1,1,1,1,2,3,4,2,3,4), ncol = 3, byrow = TRUE), widths = c(1,1, 1), heights = c(1.5, 1.5, 1.5))
  #layout(matrix(c(1, 2, 1, 3, 1, 4), ncol = 2, byrow = TRUE), widths = c(1.8, 1), heights = c(1.5, 1.5, 1.5))
  par(mar = c(5, 5, 4, 1)) 
  plot(GAM.output$`50%`~GAM.output$cov, type='l', ylim = c(0,1), lty = 1, lwd =2, col ='black', main='Generalised Additive Model', ylab = 'Occupancy', xlab = cov, 
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size)
  lines(GAM.output$`97.5%`~GAM.output$cov, lty =2)
  lines(GAM.output$`2.5%`~GAM.output$cov, lty =2)
  pos <- auto.legend.pos(y=GAM.output$`50%`, x= GAM.output$cov)
  legend(pos, legend = c('Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1,2), lwd = c(1,2), col=c('black', 'grey'), cex=1.5)
  
  linear.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Linear')
  
  mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[linear.i], "/", full.occu.mod.summary$modname[linear.i], ".rds"))$Draws$Quantiles
  
  plot(mod.output$`50%`~mod.output$cov, type='n', ylim=c(0,1), main = 'Linear Model',ylab = 'Occupancy', xlab = cov,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col='grey', border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'green3', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Linear Mean Posterior Prediction', 'GAM Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1, 1,2),, col=c('green3','black', 'grey'), cex=1.5)
  
  
  quad.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Quadratic')
  
  mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[quad.i], "/", full.occu.mod.summary$modname[quad.i], ".rds"))$Draws$Quantiles
  
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0,1), main = 'Quadratic Model',ylab = 'Occupancy', xlab = cov,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col='grey', border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'steelblue', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Quadratic Mean Posterior Prediction', 'GAM Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1, 1,2),, col=c('steelblue','black', 'grey'), cex=1.5)
  
  
  cubic.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Cubic')
  
  mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[cubic.i], "/", full.occu.mod.summary$modname[cubic.i], ".rds"))$Draws$Quantiles
  
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0,1), main = 'Cubic model',ylab = 'Occupancy', xlab = cov,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col='grey', border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'purple', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Cubic Mean Posterior Prediction', 'GAM Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1, 1,2),, col=c('purple','black', 'grey'), cex=1.5)
}

for (sp in c( 'Sus scrofa')){
  for(cov in covariates ){
    
    pdf(paste0('Outputs/Occupany_Plots', sp, '_',cov, '.pdf'), width =25, height=15)
    
    PlotAllModelsOccu(sp, cov)
    
    dev.off()
  }
}

# Abundance

PlotAllModelsAbu <- function(sp, cov){
  
  gam.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'GAM')
  
  GAM.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[gam.i], "/", full.abu.mod.summary$modname[gam.i], ".rds"))$Draws$Quantiles
  y.max <- max(GAM.output$`97.5%`) * 1.2 
  
  layout(matrix(c(1,1,1,1,1,1,2,3,4,2,3,4), ncol = 3, byrow = TRUE), widths = c(1,1, 1), heights = c(1.5, 1.5, 1.5))
  #layout(matrix(c(1, 2, 1, 3, 1, 4), ncol = 2, byrow = TRUE), widths = c(1.4, 1), heights = c(1, 1, 1))
  par(mar = c(5, 5, 4, 1)) 
  plot(GAM.output$`50%`~GAM.output$cov, type='l', ylim = c(0, y.max), lty = 1, lwd =2, col ='black', main='Generalised Additive Model', ylab = 'Abundance', xlab = cov, 
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size)
  lines(GAM.output$`97.5%`~GAM.output$cov, lty =2)
  lines(GAM.output$`2.5%`~GAM.output$cov, lty =2)
  
  pos <- auto.legend.pos(y=GAM.output$`50%`, x= GAM.output$cov)
  legend(pos, legend = c('Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1,2), lwd = c(1,2), col=c('black', 'grey'), cex=1.5)
  
  linear.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'Linear')
  
  mod.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[linear.i], "/", full.abu.mod.summary$modname[linear.i], ".rds"))$Draws$Quantiles
  
  plot(mod.output$`50%`~mod.output$cov, type='n', ylim=c(0, y.max ), main = 'Linear Model',ylab = 'Abundance', xlab = cov,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col='grey', border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'green3', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Linear Mean Posterior Prediction', 'GAM Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1, 1,2),, col=c('green3','black', 'grey'), cex=1.5)
  
  
  quad.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'Quadratic')
  
  mod.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[quad.i], "/", full.abu.mod.summary$modname[quad.i], ".rds"))$Draws$Quantiles
  
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0, y.max ), main = 'Quadratic Model',ylab = 'Abundance', xlab = cov,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col='grey', border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'steelblue', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Quadratic Mean Posterior Prediction', 'GAM Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1, 1,2),, col=c('steelblue','black', 'grey'), cex=1.5)
  
  
  cubic.i <- which(full.abu.mod.summary$Covariate == cov & full.abu.mod.summary$Species == sp & full.abu.mod.summary$Model_Type == 'Cubic')
  
  mod.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", full.abu.mod.summary$Covariate[cubic.i], "/", full.abu.mod.summary$modname[cubic.i], ".rds"))$Draws$Quantiles
  
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0, y.max ), main = 'Cubic model',ylab = 'Abundance', xlab = cov,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2,    # Increase axis label size
       cex.axis = 2,
       yaxs="i")   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col='grey', border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'purple', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Cubic Mean Posterior Prediction', 'GAM Mean Posterior Prediction', '95%-Bayesian CI'), lty = c(1, 1,2),, col=c('purple','black', 'grey'), cex=1.5)
}

for (sp in c('Macaca nemestrina', 'Sus scrofa')){
  for(cov in covariates ){
    print(sp)
    print(cov)
    pdf(paste0('Outputs/Abundance_Plots/', sp, '_',cov, '.pdf'), width =25, height=15)
    
    PlotAllModelsAbu(sp, cov)
    
    dev.off()
  }
}


# export statistics Occupancy

for (sp in c('Macaca nemestrina', 'Sus scrofa')){
  for(cov in covariates ){
    print(sp)
    print(cov)
   
   temp <- full.occu.mod.summary |>
      filter(Species == sp & Covariate == cov) |> 
      mutate_if(is.numeric, round, digits =2) |>
      select(Species,
             Covariate,
             Model_Type,
             correlation_GAM,
             nrmse_GAM,
             Bayes_P,
             C_Hat,
             DIC)
   
   write_csv(temp,
             paste0('Outputs/Model_metrics/Occupancy/', sp, '_', cov, '.csv'))
   
   temp <- full.abu.mod.summary |>
     filter(Species == sp & Covariate == cov) |>
     mutate_if(is.numeric, round, digits =2) |>
     select(Species,
            Covariate,
            Model_Type,
            correlation_GAM,
            nrmse_GAM,
            Bayes_P,
            C_Hat,
            DIC)
   
   write_csv(temp,
             paste0('Outputs/Model_metrics/Abundance/', sp, '_', cov, '.csv'))
    
  }
}









