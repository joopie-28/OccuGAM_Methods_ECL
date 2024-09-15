### Bayesian non-linear abundance model 
# version without random effect for source, geared for use in ECL dataset

#issues we have 2 lambdas now, that will not work.
# need to rename lambda to... lambda.gam?

model {
  
  # Hyperpriors for random effects
  
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
    log(lambda[i]) <- inprod(b, X[i,]) + bLandscape[Landscape[i]] + bYear[Year[i]]
    
    # Realised ecological state
    N[i] ~ dpois(lambda[i])
    
    for(j in 1:nsurvey) {
      # Detection probability modeled as a function of sampling effort i.e. camera's - note we remove source here fot ECL data
      logit(det_prob[i,j]) <- a0 + aEffort * Effort[i, j] #+ aSource[Source[i]]
      y[i,j] ~ dbin(det_prob[i,j], N[i])
      
      #### Derived Parameters
      
      # Log likelihood for WAIC 
      log_lik0[i, j] <- logdensity.bin(y[i, j], det_prob[i,j], N[i])
      
      # Bayes P-Value
      
      Presi[i,j] <- (y[i,j] - det_prob[i,j])^2       # Calculate the squared residual error of the observed data 
      y.new[i,j] ~ dbin(det_prob[i,j], N[i])              # Simulate observed data 
      Presi.new[i,j] <- (y.new[i,j] - det_prob[i,j])^2  # Calculate squared residual error of simulated data 
      
    }
    
    #### get the row log-likelihood
    log_lik[i] <- sum(log_lik0[i,])
    
  }
  
  
  SSEobs <- sum(Presi[,])     # Calculate the sum of squared residual errors for observed data
  SSEsim <- sum(Presi.new[,]) # Calculate the sum of squared residual error for the simulated data
  p.val <- step(SSEsim - SSEobs)
  
  # the detection process priors
  a0 ~ dnorm(0,0.75)
  aEffort ~ dnorm(0,0.75)
  
  ## Parametric effect priors 
  b[1] ~ dnorm(0,0.75)
  
  ## prior for s(covariate), with a maximum spline complexity of 5
  K1 <- S1[1:4,1:4] * lambda_gam[1]  + S1[1:4,5:8] * lambda_gam[2]
  b[2:5] ~ dmnorm(zero[2:5],K1) 
  
  ## smoothing parameter priors CHECK...
  for (k in 1:2) {
    lambda_gam[k] ~ dgamma(.05,.005)
    rho[k] <- log(lambda_gam[k])
  }
}
