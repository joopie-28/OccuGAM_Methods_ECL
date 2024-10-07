##### Summarise predictive accuracy of models and 

library(ubms)
library(loo)


# Bring in the models, compute WAIC etc

# Create a table

# Visualise result

# Part 1 - Load in model results

# Polynomial models

current.poly.mod <- readRDS('/Users/sassen/Desktop/Occupancy/Macaca_nemestrina_Poly_Avg_human_footprint_3km_2024-09-17_JMS.rds')

ubms::loo(current.poly.mod$Cubic)
loo::getlo(current.poly.mod$Cubic)
ubms::extract_log_lik(current.poly.mod$Cubic)

# GAM models

current.occugam <- readRDS("/Users/sassen/Desktop/02_HPC_Occupancy/results/FULL_MODELS/Macaca_nemestrina_Avg_human_footprint_3km_2024-09-16_JMS.rds")
newlist <- lapply(current.poly.mod, ubms::loo)
newlist[['GAM']] <- current.occugam$LOO
loo_compare(newlist)
extract_log_lik(current.poly.mod$Linear)

current.occugam$WAIC
loo


covs <- c("Avg_FLLI_3km","Avg_forest_cover_3km", "Avg_human_footprint_3km", "Avg_human_population_3km", "Avg_OP_percent_3km")
species <-c('Macaca nemestrina', 'Sus scrofa')



for(sp in species){
  print(sp)
  for(cov in covs){
    print(cov)
    pdf(paste0('Outputs/OccupancyModelPlots/',str_replace_all(sp," ","_"),"_",cov, 'test.pdf'), height = 10, width=20)
    #### STart visual plotting
    # Step 1 = set up plotting parameters
    par(mfrow=c(1,2))
    par(mar=c(5,5,1,1))
    
    #Terrible bit of code
    modname <- list.files(paste0("/Users/sassen/Desktop/02_HPC_Occupancy/results/PLOTTING_RDS","/", cov))[grepl(str_replace_all(sp, " ", "_"),list.files(paste0("/Users/sassen/Desktop/02_HPC_Occupancy/results/PLOTTING_RDS","/", cov)))]
    sampled.object  <- readRDS(paste0("/Users/sassen/Desktop/02_HPC_Occupancy/results/PLOTTING_RDS/",cov,"/",modname))
    full.mod<- readRDS(paste0("/Users/sassen/Desktop/02_HPC_Occupancy/results/FULL_MODELS/",modname))
    
    # Extract relevant data for easy plotting
    quants <- sampled.object$Draws$Quantiles
    draws <- sampled.object$Draws$Samples
    med.deriv <- sampled.object$Derivatives$Second$Median
    draws.deriv <- sampled.object$Derivatives$Second$Samples
    deriv.CI <- sampled.object$Derivatives$Second$Quantiles
    
    # Occupancy plotting

    # Save a better version of the variable names for plotting
    cov.name <- ifelse(cov == "Avg_human_population_3km", 'Human Population',
                       ifelse(cov == 'Avg_OP_percent_3km', 'Oil Palm Percent', 
                              ifelse(cov == "Avg_forest_cover_3km", 'Forest Cover', 
                                     ifelse(cov == "Avg_human_footprint_3km", 'Human Footprint', 'Forest Integrity'))))
    
    # Next, populate plots
    plot(`50%`~cov, data = quants,type='n', ylim=c(0,1), lwd=2, ylab = 'Occupancy', xlab = cov.name, axes=F, cex.lab = 2.5)
    
    box()
    axis(side=2, cex.axis=2)
    axis(side=1, cex.axis =2)
    
    polygon(c(quants$cov, rev(quants$cov)), 
            c(quants$`97.5%`, rev(quants$`2.5%`)), border=NA, col=alpha('lightgrey',0.2))
    
    
    lines(quants$`50%`~quants$cov, type='l', ylim=c(0,1), lwd=2.5, col='black')
    
    # Identifier plotting
    text(y=1,x=min(quants$cov),(paste0(sp, ' ~ ', cov.name)), adj=0, cex=1)
    
    gam.l.ID <- full.mod$LOO$looic
    
    legend('bottomright',
           col = 'black',
           legend = c(paste0('GAM, LOOIC = ' ,round( gam.l.ID,2))),
           lty = 2)
    
    # predict new data
    mname <- list.files('/Users/sassen/Desktop/Occupancy')[grepl(str_replace_all(paste0(sp, " Poly ", cov), " ", "_") ,list.files('/Users/sassen/Desktop/Occupancy'))]
    current.poly.mod <- readRDS(paste0('/Users/sassen/Desktop/Occupancy/',mname))
    
    
    
    model <- current.poly.mod[[1]]
    newdat <- data.frame(seq(min(umf.list.occu[[sp]]@siteCovs[cov]), max(umf.list.occu[[sp]]@siteCovs[cov]), by = 0.001))
    colnames(newdat) <- cov
  
    pred <- ubms::predict(model, submodel='state',newdata = newdat, re.form=NA)
    pred$cov <- newdat[[cov]]
    pred <- pred[order(pred$cov),]
    
    cols <- c('black','orange','steelblue')
    # Plot and make sure frame lines up with the Bayesian version
    plot(pred$Predicted~pred$cov, 
         ylab = paste0("Occupancy"),
         ylim = c(0, 1), xlab = cov.name, type='l', lwd=2.5, col=cols[1],
         cex.lab = 2.5, cex.axis =2)
    polygon(c(pred$cov, rev(pred$cov)),
            c(pred$`97.5%`, rev(pred$`2.5%`)),  border=NA, col=alpha(cols[i],0.2))
    
    for (i in c(2,3)){
      model <- current.poly.mod[[i]]
      newdat <- data.frame(seq(min(umf.list.occu[[sp]]@siteCovs[cov]), max(umf.list.occu[[sp]]@siteCovs[cov]), by = 0.001))
      colnames(newdat) <- cov
      
      pred <- ubms::predict(model, submodel='state',newdata = newdat, re.form=NA)
      pred$cov <- newdat[[cov]]
      pred <- pred[order(pred$cov),]
      lines(pred$Predicted~pred$cov,lwd=2.5, col=cols[i])
      # Add the intervals
      polygon(c(pred$cov, rev(pred$cov)),
              c(pred$`97.5%`, rev(pred$`2.5%`)),  border=NA, col=alpha(cols[i],0.2))
    }
    
    # Fit estimates
    l.IC <-lapply(current.poly.mod, ubms::loo)
    legend('bottomright',
           col = cols,
           legend = c(paste0('Linear, LOOIC = ' ,round(l.IC[[1]]$looic,2)),
                      paste0('Quadratic, LOOIC = ' ,round(l.IC[[2]]$looic,2)),
                     paste0('Cubic, LOOIC = ' ,round(l.IC[[3]]$looic,2))),
           lty = 2)
  
    dev.off()
    # Add the intervals
  }
}

