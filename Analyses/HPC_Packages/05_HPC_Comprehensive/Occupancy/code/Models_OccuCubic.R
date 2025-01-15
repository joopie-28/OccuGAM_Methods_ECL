# Bayesian Occupancy model with cubic function of covariates (third degree polynomial)

# we will need to create new variables for model fitting

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
    # Latent state, b is the intercept + smoothing term, bLandscape is the random effect on the intercept and bYear is the year random effect on the intercept
    logit(psi[i]) <- inprod(b, DisCov[i,]) + bLandscape[Landscape[i]] + bYear[Year[i]]
    
    # Realised ecological state
    z[i] ~ dbern(psi[i])
    
    for(j in 1:nsurvey) {
      # Detection probability modelled as a function of sampling effort i.e. camera's
      logit(det_prob[i,j]) <- lp[i,j]
      mu.lp[i, j] <- a0 + aEffort * Effort[i, j]
      
      lp[i, j] ~ dnorm(mu.lp[i, j], tau.p)
      y[i,j] ~ dbern(det_prob[i,j] * z[i])
      
      #### Derived Parameters
      # Bayes P-Value
      
      Presi[i,j] <- (y[i,j] - det_prob[i,j])^2       # Calculate the squared residual error of the observed data 
      y.new[i,j] ~ dbern(det_prob[i,j] * z[i])              # Simulate observed data 
      Presi.new[i,j] <- (y.new[i,j] - det_prob[i,j])^2  # Calculate squared residual error of simulated data 
      
      # log-likelihood for WAIC (Occupancy)
      log_lik0[i,j] <- logdensity.bern(y[i,j], det_prob[i,j] * z[i])
    }
    
    # Calculate row-level (site level) log-likelihood
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
  b[1] ~ dnorm(0,0.75) 
  b[2] ~ dnorm(0,0.75) 
  b[3] ~ dnorm(0,1)
  b[4] ~ dnorm(0,1) 
  tau.p <- pow(sd.p, -2)
  sd.p ~ dunif(0, 3)
  
}
