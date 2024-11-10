##########################
### Analyse HPC output ###
##########################

# GAMs for Wildlife
# 10/11/2024
# J.M. Sassen

# Packages
library(data.table)
library(hydroGOF)

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







#### Plotting #####
plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0,100))
lines(mod.output$`97.5%`~mod.output$cov)
lines(mod.output$`2.5%`~mod.output$cov)

plot(GAM.output$`50%`~GAM.output$cov, type='l', ylim = c(0,100))






lines(mod.extract$Draws$Quantiles$`50%`~ mod.extract$Draws$Quantiles$cov, type='l')




# SOMETHING MESSY WITH THE FITTING PROCESS! - need to fix now.

test.pcount <- unmarked::pcount(~num_cams_active_at_date~poly(Avg_FLLI_3km,3) + (1|Landscape) , data = umf.list.abu$`Macaca nemestrina`)


# Make predictions based on the fitted model
predictions <- unmarked::predict(test.pcount, type = "state", re.form=NA )

# Combine predictions with the footprint values for plotting
plot_data <- data.frame(Avg_FLLI_3km = umf.list.abu$`Macaca nemestrina`@siteCovs$Avg_FLLI_3km, 
                        predicted_abundance = predictions$Predicted)

plot(plot_data$predicted_abundance~plot_data$Avg_FLLI_3km, type='p', ylim=c(0,200))
# Plot the results
ggplot(plot_data, aes(x = Avg_human_footprint_3km, y = predicted_abundance)) +
  geom_line() +
  labs(x = "Average Human Footprint (3km)", 
       y = "Predicted Abundance of Sus scrofa", 
       title = "Predicted Abundance of Sus scrofa vs. Human Footprint") +
  theme_minimal() +ylim(0,35)


GAM.output <- readRDS(paste0(result.folder, "/Abundance/results/PLOTTING_RDS/", cov, "/", GAM.mod, ".rds"))$Draws$Quantiles

lines(GAM.output$`50%`~GAM.output$cov )
lines(GAM.output$`2.5%`~GAM.output$cov )
lines(GAM.output$`97.5%`~GAM.output$cov )
tmp_est <- exp(intercept + as.matrix(newdat[c('DisCov', 'DisCov2', 'DisCov3')]) %*% t(mod_mcmc[,2:4]))

lines(mod.extract$Draws$Quantiles$`50%` ~ mod.extract$Draws$Quantiles$cov, col = 'red', lwd =1.5) 
points(mod$Full_model$q50$lambda~mod$Data_List$DisCov)


toplot <- matrix(NA, 1000, 3)
for(i in 1:1000) {
  p.tmp <- with(mod$Full_model$sims.list,
                      b0 + b1 * newdat[i,'DisCov'] + b2 * newdat[i,'DisCov2'] + b3 * newdat[i,'DisCov3'])
  p.tmp <- exp(p.tmp)
  toplot[i, 1] <- mean(p.tmp)
  toplot[i, 2:3] <- hdi(p.tmp)
}

plot(toplot[,1]~newdat$DisCov, ylim=range(toplot), type='l', las=1,
     xlab="HFP", ylab="ABU")

polygon(x=c(xx, rev(xx)),
        y=c(toplot[, 2], rev(toplot[, 3])),
        col=adjustcolor('skyblue', 0.5), border=NA)

mod$Full_model$q50


plot(mod.extract$Draws$Quantiles$`50%`~ mod.extract$Draws$Quantiles$cov)

plot(mod$Full_model$q50$N~ mod$Data_List$DisCov, col=as.factor(mod$Data_List$Landscape))


exp(mod$Full_model$q50$b0 * rep(1, 524) + 
mod$Full_model$q50$b1 * mod$Data_List$DisCov +
mod$Full_model$q50$b2 * mod$Data_List$DisCov2 +
mod$Full_model$q50$b3 * mod$Data_List$DisCov3 +
  mod$Full_model$q50$bLandscape)



