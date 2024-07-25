# Function for plotting occupancy models fitted in JAGS and Unmarked
GenerateModelPlots <- function(sampled.object, plot.draws = F, plot.freq.occu = F, plot.deriv = F, 
                                   cov.name, additional.plot = NULL,fit.col,offset =NULL){
  
  cov.name <- ifelse(cov.name == "Avg_human_population_3km", 'Human Population',
                     ifelse(cov.name == 'Avg_OP_percent_3km', 'Oil Palm Percent', 
                            ifelse(cov.name == "Avg_forest_cover_3km", 'Forest Cover', 
                                   ifelse( cov.name == "Avg_human_footprint_3km", 'Human Footprint', 'Forest Integrity'))))
  
  quants <- sampled.object$Draws$Quantiles
  draws <- sampled.object$Draws$Samples
  
  med.deriv <- sampled.object$Derivatives$Second$Median
  draws.deriv <- sampled.object$Derivatives$Second$Samples
  
  deriv.CI <- sampled.object$Derivatives$Second$Quantiles
  
  # Step 1 = set up plotting parameters
  
  
  par( mar=c(5,5,1,5))
  # Next, populate plots
  plot(`50%`~cov, data = quants,type='n', ylim=c(0,1), lwd=2, ylab = 'Occupancy', xlab = cov.name, axes=F, cex.lab = 2.5)
  
  box()
  axis(side=2, cex.axis=2)
  axis(side=1, cex.axis =2)
  
  polygon(c(quants$cov, rev(quants$cov)), 
          c(quants$`94.5%`, rev(quants$`5.5%`)), border=NA, col=alpha('lightgrey',0.2))
  
  
  lines(quants$`50%`~quants$cov, type='l', ylim=c(0,1), lwd=2.5, col=fit.col)
  
  # Identifier plotting
  edf <- sampled.object$EDF
  species <- sampled.object$Species
  
  
  text(y=.85,x=min(quants$cov)+1 + offset, paste0('EDF = ~ ', round(edf,2)), adj=0, cex=2.5)
  text(y=.95,x=min(quants$cov)+1 + offset, substitute(italic(species)), adj=0, cex=2.5)
  
  
  plotting.list <- additional.plot$plotting
  
  # add trasngression
  for(i in 1:length(plotting.list)){
    
    if(!is.null(plotting.list[[i]])){
      
      lines(plotting.list[[i]]$seq.fitted ~ plotting.list[[i]]$seq.cov, col=fit.col, lwd=8)
      
    }
  }
  
  
  
}
