
# Think we should show all models
PlotBayesianComparisonModels <- function(sp, cov, full.mod.summary,
                                 rds.path, umf.list.abu, 
                                 umf.list.occu){
  # Extract relevant information
  #  sp <- "Macaca nemestrina"
  #  cov <- "Avg_OP_percent_3km"
  
  # Subset the relevant rows from our total dataframe so we have the modname we used
  # to store the JAGS model file.
  current.modname <- full.mod.summary$ModName[full.mod.summary$Species == sp & full.mod.summary$Covariate == cov]
  
  # Isolate the plotting object from the results folder
  sampled.object <- readRDS(paste0(rds.path, "/PLOTTING_RDS/", cov, "/", current.modname)) 
  
  #GenerateModelPlots(current.plot, plot.draws = F, type = 'abundance', 
  #                   cov.name = cov, fit.col = 'black')
  
  # Set up meta-plotting environment to accomodate 2 screens
  par(mfrow = c(1,3))
  
  # Plotting JAGS Models
  
  # Save a better version of the variable names for plotting
  cov.name <- ifelse(cov == "Avg_human_population_3km", 'Human Population',
                     ifelse(cov == 'Avg_OP_percent_3km', 'Oil Palm Percent', 
                            ifelse(cov == "Avg_forest_cover_3km", 'Forest Cover', 
                                   ifelse(cov == "Avg_human_footprint_3km", 'Human Footprint', 'Forest Integrity'))))
  
  # Extract relevant data for easy plotting
  quants <- sampled.object$Draws$Quantiles
  draws <- sampled.object$Draws$Samples
  med.deriv <- sampled.object$Derivatives$Second$Median
  draws.deriv <- sampled.object$Derivatives$Second$Samples
  deriv.CI <- sampled.object$Derivatives$Second$Quantiles
  
  # Occupancy plotting
  if(type == 'occupancy'){
    
    # Step 1 = set up plotting parameters
    
    par(mar=c(5,5,1,1))
    # Next, populate plots
    plot(`50%`~cov, data = quants,type='n', ylim=c(0,1), lwd=2, ylab = 'Occupancy', xlab = cov.name, axes=F, cex.lab = 2.5)
    
    box()
    axis(side=2, cex.axis=2)
    axis(side=1, cex.axis =2)
    
    polygon(c(quants$cov, rev(quants$cov)), 
            c(quants$`97.5%`, rev(quants$`2.5%`)), border=NA, col=alpha('lightgrey',0.2))
    
    
    lines(quants$`50%`~quants$cov, type='l', ylim=c(0,1), lwd=2.5, col='black')
    
    # Identifier plotting
    text(y=.95,x=min(quants$cov),(paste0(sp, ' ~ ', cov.name)), adj=0, cex=1)
    
  }
  
  # Abundance plotting
  if(type == 'abundance'){
    
    # Find the maximum abundance so we can set up a reasonable window
    max.abu <- max(sampled.object$Draws$Quantiles$`97.5`)
    
    # Step 1 = set up plotting parameters
    
    par(mar=c(5,5,1,1))
    # Next, populate plots
    plot(`50%`~cov, data = quants,type='n',ylim=c(0, (max.abu)), lwd=2, ylab = 'Abundance', xlab = cov.name, axes=F, cex.lab = 2.5)
    
    box()
    axis(side=2, cex.axis=2)
    axis(side=1, cex.axis =2)
    
    polygon(c(quants$cov, rev(quants$cov)), 
            c(quants$`97.5%`, rev(quants$`2.5%`)), border=NA, col=alpha('lightgrey',0.2))
    
    
    lines(quants$`50%`~quants$cov, type='l', ylim=c(0,1), lwd=2.5, col='black')
    
    # Identifier plotting
    
    text(y=max.abu-0.5,x=min(quants$cov),(paste0(sp, ' ~ ', cov.name)), adj=0, cex=2.5)
    
  }
  
  # Fit and plot the frequentists versions using unmarked with increasing polynomial terms of the covariate
  models <- list()
  formulae <- list(
    "Linear" = as.formula(paste("~ num_cams_active_at_date ~", cov,"+ (1|Landscape) + (1|Year)")),
    "Quadratic" = as.formula(paste("~ num_cams_active_at_date ~ poly(", cov, ", 2, raw = TRUE) + (1|Landscape) + (1|Year)")),
    "Cubic" = as.formula(paste("~ num_cams_active_at_date ~ poly(", cov, ", 3, raw = TRUE) + (1|Landscape) + (1|Year)")))
  
  
  # Fit the models
  
  # Occupancy pathway
  if(type == 'occupancy'){
    
    # Set up the year covariate for the data
    
    umf.list.occu[[sp]]@siteCovs$Year <- as.numeric(as.factor(stringr::str_extract(umf.list.occu[[sp]]@siteCovs[['survey_id']], 
                                                                                   stringr::regex("(\\d+)(?!.*\\d)"))))
    
    # occu models
    for (name in names(formulae)) {
      print(paste0('Fitting frequentist models, currently on the ', name, ' model'))
      
      #occu models
      models[[name]] <- ubms::stan_occu(formulae[[name]], data = umf.list.occu[[sp]])
    }
    
    
    # Select the best model
    best.model <- models[[aictab(models)[1,1]]]
    degree <- aictab(models)[1,1]
    
    # Plot the best model next to the Bayesian version
    
    pred <- ubms::predict(best.model, type='state', re.form=NA)
    pred$cov <- umf.list.occu[[sp]]@siteCovs[[cov]]
    pred <- pred[order(pred$cov),]
    
    # Plot and make sure frame lines up with the Bayesian version
    plot(pred$Predicted~pred$cov, 
         ylab = paste0("Occupancy (",degree,")" ),
         ylim = c(0, 1), xlab = cov.name, type='l', lwd=2.5, col='black',
         cex.lab = 2.5, cex.axis =2 )
    
    # Add the intervals
    polygon(c(pred$cov, rev(pred$cov)),
            c(pred$upper, rev(pred$lower)),  border=NA, col=alpha('lightgrey',0.2))
    
    # The final panel is the difference graph, which shows the subtrated values + sample size!
    
    # Subtract predictions
    diff <- (quants$`50%` - pred$Predicted)
    plot(diff ~ quants$cov,  ylim = c(min(diff), c(max(diff))), xlab = cov.name, type='l', lwd=2.5, col='black',
         cex.lab = 2.5, cex.axis =2, lty=2 , ylab = 'Difference between Polynomial and Bayesian GAM')
    
    # Add bars to show the dataspread
    
    samples.mat <- as.matrix(table(quants$cov))
    
    # Sample data: x-values, y-values, and sample counts
    x_values <- as.numeric(rownames(samples.mat))
    sample_counts <- rep(1, length(samples.mat[,1])) # Number of samples taken at each x position
    
    # Add bars at each x-axis value to indicate sample counts
    segments(x0 = x_values, y0 = 0, x1 = x_values, y1 = sample_counts, 
             col = alpha("black",0.2), lwd = 3)
    
  }
  
  # Abundance pathway
  if(type == 'abundance'){
    
    # Set up the year covariate for the data
    
    umf.list.abu[[sp]]@siteCovs$Year <- as.numeric(as.factor(stringr::str_extract(umf.list.abu[[sp]]@siteCovs[['survey_id']], 
                                                                                  stringr::regex("(\\d+)(?!.*\\d)"))))
    
    # Abundance models
    for (name in names(formulae)) {
      print(paste0('Fitting frequentist models, currently on the ', name, ' model'))
      
      # pcount models = RN
      models[[name]] <- tryCatch({
        ubms::stan_pcount(formulae[[name]], data = umf.list.abu[[sp]])
      }, error = function(e) {
        message(paste("Error in fitting model", name, ":", e$message))
        return(NULL)  # Return NULL when an error occurs
      })
    }
    
    
    # Select the best model
    best.model <- models[[aictab(models)[1,1]]]
    degree <- aictab(models)[1,1]
    
    # Plot the best model next to the Bayesian version
    
    pred <- ubms::predict(best.model, type='state', re.form=NA)
    pred$cov <- umf.list.abu[[sp]]@siteCovs[[cov]]
    pred <- pred[order(pred$cov),]
    
    # Plot and make sure frame lines up with the Bayesian version
    plot(pred$Predicted~pred$cov, 
         ylab = paste0("Abundance (",degree,")" ),
         ylim = c(0, max.abu), xlab = cov.name, type='l', lwd=2.5, col='black',
         cex.lab = 2.5, cex.axis =2 )
    
    # Add the intervals
    polygon(c(pred$cov, rev(pred$cov)),
            c(pred$upper, rev(pred$lower)),  border=NA, col=alpha('lightgrey',0.2))
    
    # The final panel is the difference graph, which shows the subtrated values + sample size!
    
    # Subtract predictions
    diff <- (quants$`50%` - pred$Predicted)
    
    # Clever way to set y axis for the difference
    small.diff <- min(diff)
    large.diff <- max(diff)
    
    if(abs(large.diff) < 2 & abs(small.diff) < 2){
      small.diff <- -2
      large.diff <- 2
    }
    
    
    plot(diff ~ quants$cov,  ylim = c(small.diff , large.diff), xlab = cov.name, type='l', lwd=2.5, col='black',
         cex.lab = 2.5, cex.axis =2, lty=2 , ylab = 'Difference between Polynomial and Bayesian GAM')
    
    # Add bars to show the dataspread
    
    samples.mat <- as.matrix(table(quants$cov))
    
    # Sample data: x-values, y-values, and sample counts
    x_values <- as.numeric(rownames(samples.mat))
    sample_counts <- rep(1, length(samples.mat[,1])) # Number of samples taken at each x position
    
    # Add bars at each x-axis value to indicate sample counts
    points(x = x_values, y = rep(0, length(x_values)), col = alpha("black",0.2))
    
  }
  
}
