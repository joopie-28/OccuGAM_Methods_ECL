PlotAllModelsOccu <- function(sp, cov){
  cov.name <- ifelse(cov == "Avg_human_population_3km", 'Human Population',
                     ifelse(cov == 'Avg_OP_percent_3km', 'Oil Palm Percent', 
                            ifelse(cov == "Avg_forest_cover_3km", 'Forest Cover (%)', 
                                   ifelse(cov == "Avg_human_footprint_3km", 'Human Footprint (%)', 'Forest Integrity'))))
  
  # revert standardisation for plotting
  mean <- sp_standardisation_params[[sp]][sp_standardisation_params[[sp]]$Variable == cov,2]
  sd <- sp_standardisation_params[[sp]][sp_standardisation_params[[sp]]$Variable == cov,3]
  
  gam.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'GAM')
  
  GAM.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[gam.i], "/", full.occu.mod.summary$modname[gam.i], ".rds"))$Draws$Quantiles
  
  # back transform
  GAM.output$cov <- ((GAM.output$cov * sd) + mean)
  y.max <- 1
  
  layout(matrix(c(1,2,3), ncol = 3, byrow = TRUE), widths = c(1,1, 1,1), heights = c(1.5, 1.5, 1.5, 2.25,2.25))
  par(mar = c(5, 5.5, 4, 1)) 
  
  linear.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Linear')
  
  mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[linear.i], "/", full.occu.mod.summary$modname[linear.i], ".rds"))$Draws$Quantiles
  # back transform
  mod.output$cov <- ((mod.output$cov * sd) + mean)
  
  plot(mod.output$`50%`~mod.output$cov, type='n', ylim=c(0, y.max), main = 'Linear versus GAM',ylab = '', xlab = cov.name,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2.5,    # Increase axis label size
       cex.axis = 2.5,
       yaxs="i",
       las=1,
       mgp = c(3.65,1.4,0))   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col=alpha('grey',0.4), border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'green3', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Linear', 'GAM', '95%-CI'), lty = c(1, 1,2), col=c('green3','black', 'grey'), cex=1.5)
  
  
  quad.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Quadratic')
  
  mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[quad.i], "/", full.occu.mod.summary$modname[quad.i], ".rds"))$Draws$Quantiles
  # back transform
  mod.output$cov <- ((mod.output$cov * sd) + mean)
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0, y.max ), main = 'Quadratic versus GAM',ylab = '',xlab = cov.name,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2.5,    # Increase axis label size
       cex.axis = 2.5,
       yaxs="i",
       las=1,
       mgp = c(3.65,1.4,0))   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col=alpha('grey',0.4), border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'steelblue', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Quadratic', 'GAM', '95%-CI'), lty = c(1, 1,2),, col=c('steelblue','black', 'grey'), cex=1.5)
  
  
  cubic.i <- which(full.occu.mod.summary$Covariate == cov & full.occu.mod.summary$Species == sp & full.occu.mod.summary$Model_Type == 'Cubic')
  
  mod.output <- readRDS(paste0(result.folder, "/Occupancy/results/PLOTTING_RDS/", full.occu.mod.summary$Covariate[cubic.i], "/", full.occu.mod.summary$modname[cubic.i], ".rds"))$Draws$Quantiles
  # back transform
  mod.output$cov <- ((mod.output$cov * sd) + mean)
  
  plot(mod.output$`50%`~mod.output$cov, type='l', ylim=c(0, y.max ), main = 'Cubic versus GAM',ylab = '', xlab = cov.name,
       cex.main = 2.5,   # Increase title size
       cex.lab = 2.5,    # Increase axis label size
       cex.axis = 2.5,
       yaxs="i",
       las=1,
       mgp = c(3.65,1.4,0))   # Increase axis tick label size))
  lines(mod.output$`97.5%`~mod.output$cov, lty =2)
  lines(mod.output$`2.5%`~mod.output$cov, lty =2)
  polygon(c(mod.output$cov, rev(mod.output$cov)), c(mod.output$`50%`, rev(GAM.output$`50%`)), col=alpha('grey',0.4), border =NA)
  lines(GAM.output$`50%`~GAM.output$cov,lty = 1, lwd =2, col ='black')
  lines(mod.output$`50%`~mod.output$cov, lty =1,col = 'purple', lwd=2)
  pos <- auto.legend.pos(y=mod.output$`50%`, x= mod.output$cov)
  legend(pos, legend = c('Cubic', 'GAM', '95%-CI'), lty = c(1, 1,2),, col=c('purple','black', 'grey'), cex=1.5)
}