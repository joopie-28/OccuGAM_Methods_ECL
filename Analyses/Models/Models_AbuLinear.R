# Bayesian N-Mixture model with linear function of covariates (first degree polynomial)


model {
  
  # There are 2 random effects, we define those here
  
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
    # Latent state (abundance), b is the intercept + smoothing term, bLandscape is the random effect on the intercept and bYear is the year random effect on the intercept
    log(lambda[i]) <-  b0 + b1 * DisCov[i] + bLandscape[Landscape[i]] + bYear[Year[i]]
    
    # Realised ecological state
    N[i] ~ dpois(lambda[i])
    
    
    for(j in 1:nsurvey) {
      # Detection probability modeled as a function of sampling effort i.e. camera's - note we remove source here fot ECL data
      logit(det_prob[i,j]) <- a0 + aEffort * Effort[i, j] #+ aSource[Source[i]]
      y[i,j] ~ dbin(det_prob[i,j], N[i])
      
      #### Derived Parameters
      # Bayes P-Value
      
      Presi[i,j] <- (y[i,j] - det_prob[i,j])^2       # Calculate the squared residual error of the observed data 
      y.new[i,j] ~ dbern(det_prob[i,j] * z[i])              # Simulate observed data 
      Presi.new[i,j] <- (y.new[i,j] - det_prob[i,j])^2  # Calculate squared residual error of simulated data 
      
      # Log likelihood for WAIC 
      log_lik0[i, j] <- logdensity.bin(y[i, j], det_prob[i,j], N[i])
    }
    #### get the row log-likelihood
    log_lik[i] <- sum(log_lik0[i,])
  }
  
  SSEobs <- sum(Presi[,])     # Calculate the sum of squared residual errors for observed data
  SSEsim <- sum(Presi.new[,]) # Calculate the sum of squared residual error for the similuated data
  p.val <- step(SSEsim - SSEobs)
  c.hat <- SSEsim / SSEobs
  
  # the detection process priors
  a0 ~ dnorm(0,0.75)
  aEffort ~ dnorm(0,0.75)
  
  ## Parametric effect priors 
  
  b0 ~ dunif(-5, 5) 
  b1 ~ dunif(-5, 5) 
  
  
}
