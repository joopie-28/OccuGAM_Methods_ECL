#######################################################################
### Non-Linear responses to environmental covariates in SEA Mammals ###
#######################################################################

# Methods comparison paper

# J.M. Sassen 
# Ecological Cascades Lab, University of Queensland
# 12-07-2024

# Main Analysis Part 3: Inspect and visualise model results

#### Activate Packages


#### Source Functions



#### Produce the plots

# take the differences of all species nonlinear wrt to their linears

#### Import the results
rds.path <- "/Users/sassen/Dropbox/Joop MSc QBIO project/Wildlife GAMs Methods Paper/GAMS Methods - Results/Abundance Models/01_HPC_Abundance_v300724/results" # should be a dropbox folder

# remember that we have subfolders for each variable
covariates <-  c("Avg_FLLI_3km", 
                 "Avg_forest_cover_3km", 
                 "Avg_human_footprint_3km", 
                 "Avg_OP_percent_3km", 
                 "Avg_human_population_3km")

# Read in all the model artefacts
mod.summary.list <- list()

# Loop through the folders to extract our model info, its broken up 
# into different 1-row df's so we need to extract and bind.
for(cov in covariates){
  current.folder <- list.files(paste0(rds.path, "/ABUGAM_RDS/", cov)) 
  
  for(mod in current.folder){
    print(mod)
    temp <- readRDS(paste0(rds.path, "/ABUGAM_RDS/", cov, "/", mod)) 
    slurm.index <- which(temp$Species != "")
    mod.summary.list[[slurm.index]] <- temp[slurm.index,]
  }
}

# Final Modelling Outputs
full.mod.summary <- rbindlist(mod.summary.list)

# Clean up environment
rm(cov, mod, mod.summary.list,slurm.index, temp)

# Produce non-linearity plots

full.mod.summary$ModName <- NA

for(cov in covariates){
  
  current.folder <- list.files(paste0(rds.path, "/PLOTTING_RDS/", cov)) 
  
  for(mod in current.folder){
    
    #GenerateThresholdPlots(current.plot, plot.draws = F, plot.freq.occu = F, plot.deriv=T, cov.name = cov)
    
    # First check if the derivative was actually approximated
    current.plot <- readRDS(paste0(rds.path, "/PLOTTING_RDS/", cov, "/", mod)) 
    
    # perform binning 
    sp <- current.plot$Species
    row.ind <- which(full.mod.summary$Species == sp & full.mod.summary$Covariate == cov)
    full.mod.summary$ModName[row.ind] <- mod
    rm(current.plot)
    
  }
}

# plot
for (i in 1:nrow(full.mod.summary)){
  mod <- full.mod.summary$ModName[i] 
  cov <-full.mod.summary$Covariate[i] 
  
  current.plot <- readRDS(paste0(rds.path, "/PLOTTING_RDS/", cov, "/", mod)) 
  GenerateModelPlots(current.plot, plot.draws = F, type = 'abundance', 
                     cov.name = cov, fit.col = 'black')
  Sys.sleep(2)
}

nrow(full.mod.summary)
#### Fit frequentist linear models
Frequentist.Model.List <- list()
for(i in c(8,10:40)){
  print(i)
  # extract species
  sp <- full.mod.summary$Species[i]
  current.cov <- full.mod.summary$Covariate[i]
  UMF.data <- umf.list.abu[[sp]]
  
  # fit the models
  Frequentist.Model.List[[paste0(sp,'-',current.cov)]] <- fitPcountMods(UMF.data, cov = current.cov)
  
}

saveRDS(Frequentist.Model.List, "Outputs/FreqModels.rds")
## fill in
Frequentist.Model.List <-append(Frequentist.Model.List, "MN_FC", after = 6)
Frequentist.Model.List <-append(Frequentist.Model.List, "MN_HP", after = 8)

# plot next to each other

pdf("Outputs/Abundance_Model_Plots/FullModelResults.pdf",
    width = 15,
    height=10)

# Set columns
par(mfrow=c(1,2), mar=c(1,1,1,1))

for (i in 1:6){
  mod <- full.mod.summary$ModName[i] 
  cov <-full.mod.summary$Covariate[i] 
  
  current.plot <- readRDS(paste0(rds.path, "/PLOTTING_RDS/", cov, "/", mod)) 
  GenerateModelPlots(current.plot, plot.draws = F, type = 'abundance', 
                     cov.name = cov, fit.col = 'black')
  # en hier de freq
  current.pred <- Frequentist.Model.List[[i]]$Prediction
  plotmax <- max(current.pred$upper) +2
  plot(current.pred$Predicted~current.pred$cov, ylab = "Abundance (Linear)",ylim=c(0,plotmax), xlab = cov, type='l')
  polygon(c(current.pred$cov, rev(current.pred$cov)),
          c(current.pred$upper, rev(current.pred$lower)),  border=NA, col=alpha('lightgrey',0.2))
  sp<- full.mod.summary$Species[i]
  
  legend("topleft", legend = sp, cex = 2.5, bty = "n")

}

for (i in c(7,9)){
  mod <- full.mod.summary$ModName[i] 
  cov <-full.mod.summary$Covariate[i] 
  
  current.plot <- readRDS(paste0(rds.path, "/PLOTTING_RDS/", cov, "/", mod)) 
  GenerateModelPlots(current.plot, plot.draws = F, type = 'abundance', 
                     cov.name = cov, fit.col = 'black')
  # en hier de freq
 
  
  plot(current.pred$Predicted~current.pred$cov, ylab = "Abundance (Linear)",ylim=c(0,50), xlab = cov, type='n')
 
  sp<- full.mod.summary$Species[i]
  
  legend("topleft", legend = sp, cex = 2.5, bty = "n")
  
}



for (i in 10:40){
  mod <- full.mod.summary$ModName[i] 
  cov <-full.mod.summary$Covariate[i] 
  
  current.plot <- readRDS(paste0(rds.path, "/PLOTTING_RDS/", cov, "/", mod)) 
  GenerateModelPlots(current.plot, plot.draws = F, type = 'abundance', 
                     cov.name = cov, fit.col = 'black')
  # en hier de freq
  current.pred <- Frequentist.Model.List[[i]]$Prediction
  plotmax <- max(current.pred$upper) +2
  plot(current.pred$Predicted~current.pred$cov, ylab = "Abundance (Linear)",ylim=c(0,plotmax), xlab = cov, type='l')
  polygon(c(current.pred$cov, rev(current.pred$cov)),
          c(current.pred$upper, rev(current.pred$lower)),  border=NA, col=alpha('lightgrey',0.2))
  sp<- full.mod.summary$Species[i]
  
  legend("topleft", legend = sp, cex = 2.5, bty = "n")
}

dev.off()



#########################

# Use function 9 to fit the UBMS model and choose the right ones.
# JAGS Models come from the HPC and are imported using function 9
# Still need to compute AUC or something, but at least we gave GOF using LOO 

# import a mod 

temp.mod <- readRDS("/Users/sassen/Desktop/01_HPC_Abundance/results/FULL_MODELS/Macaca_nemestrina_Avg_FLLI_3km_2024-09-12_JMS.rds")
temp.mod$WAIC
temp.mod$LOO

models$Linear@loo


test@submodels@submodels$state@type
summary(test)




