### Bayesian non-linear abundance model 
# version without random effect for source, geared for use in ECL dataset

# Final Bayesian Occupancy GAMM with randommeffects tructure for landscape, year and survey source

model {
  
  # There are 3 random effects, we define those here
  
  # Landscape RE hyper priors 
  sigma.bLandscape ~ dunif(0,1)
  var.bLandscape <- 1/(sigma.bLandscape*sigma.bLandscape) 
  
  # Set up random effects VARIATION FROM MEAN per Landscape
  for(k in 1:nLandscape){
    bLandscape[k] ~ dnorm(0,var.bLandscape)
  }
  
  # Year RE hyper priors 
  sigma.bYear ~ dunif(0,1)
  var.bYear <- 1/(sigma.bYear*sigma.bYear) 
  
  # Set up random effects VARIATION FROM MEAN per year
  for(k in 1:nYear){
    bYear[k] ~ dnorm(0,var.bYear)
  }
  
  
  for(i in 1:nsite) {
    # Latent state, b is the intercept + smoothing term, bLandscape is the random effect on the intercept and bYear is the year random effect on the intercept
    logit(psi[i]) <- inprod(b,X[i,]) + bLandscape[Landscape[i]] + bYear[Year[i]]
    
    # Realised ecological state
    z[i] ~ dbern(psi[i])
    
    for(j in 1:nsurvey) {
      # Detection probability modelled as a function of sampling effort i.e. camera's
      logit(det_prob[i,j]) <- a0 + aEffort * Effort[i, j] 
      y[i,j] ~ dbern(det_prob[i,j] * z[i])
      
      # log-likelihood for WAIC (Occupancy)
      log_lik0[i,j] <- logdensity.bern(y[i,j], det_prob[i,j] * z[i])
      
      #### Derived Parameters
      #Bayes P-Value
      
      Presi[i,j] <- (y[i,j] - det_prob[i,j])^2       # Calculate the squared residual error of the observed data 
      y.new[i,j] ~ dbern(det_prob[i,j] * z[i])              # Simulate observed data 
      Presi.new[i,j] <- (y.new[i,j] - det_prob[i,j])^2  # Calculate squared residual error of simulated data 
      
    }
    
    # Calculate row-level (site level) log-likelihood
    log_lik[i] <- sum(log_lik0[i,])
    
  }
  
  SSEobs <- sum(Presi[,])     # Calculate the sum of squared residual errors for observed data
  SSEsim <- sum(Presi.new[,]) # Calculate the sum of squared residual error for the similuated data
  p.val <- step(SSEsim - SSEobs)
  
  # the detection process priors
  a0 ~ dnorm(0,0.75)
  aEffort ~ dnorm(0,0.75)
  ## Parametric effect priors 
  b[1] ~ dnorm(0,0.75)
  ## prior for s(covariate), with a maximum spline complexity of 5
  K1 <- S1[1:4,1:4] * lambda[1]  + S1[1:4,5:8] * lambda[2]
  b[2:5] ~ dmnorm(zero[2:5],K1) 
  ## smoothing parameter priors CHECK...
  for (k in 1:2) {
    lambda[k] ~ dgamma(.05,.005)
    rho[k] <- log(lambda[k])
  }
}
