#### Script to run 9,600 MVGAM models in parallel in a HPC environment

# J.M. Sassen 08-03-2025

# Install packages
library(mvgam)
library(dplyr)
library(patchwork)
library(scoringRules)

##### Read in functions here
sim_piecewise <- function(x) {
  out <- vector(length = length(x))
  q_0.35 <- quantile(x, probs = 0.35)
  for (i in seq_along(x)) {
    if (x[i] <= q_0.35) {
      out[i] <- 0.75 * x[i]
    } else {
      out[i] <- (0.025 * x[i]) +
        (0.725 * q_0.35)
    }
  }
  return(as.vector(scale(out)))
}
sim_monotonic <- function(x) {
  out <- exp(2.2 * x) / (6 + exp(2 * x)) * -1
  return(as.vector(scale(out)))
}
sim_nonmonotonic <- function(x) {
  as.vector(scale(
    -0.00005 * x^11 -
      0.00020 * x^10 -
      0.00010 * x^9  +
      0.00080 * x^8  -
      0.00600 * x^7  -
      0.05000 * x^6  +
      0.07000 * x^5  +
      0.40000 * x^4  -
      0.30000 * x^3  -
      1.30000 * x^2  +
      1.3000 * x
  ))
}
sim_cubic <- function(x) {
  out <- -0.5 * x^2 + 1.5 * x^3 * -1
  return(as.vector(scale(out)))
}
sim_linear <- function(x) {
  out <- 0.4 * x
  return(as.vector(scale(out)))
}
sim_poly_then_exp <- function(x) {
  x_thresh <- .1
  
  # Centered covariate for smooth polynomial dip
  poly_part <- 10 * (x - x_thresh)^4  # Quartic, starts at ~0.5, dips
  
  # Exponential increase after threshold
  exp_part <- exp(4 * (x - x_thresh)) 
  exp_vals <- exp_part[x > x_thresh]
  
  # Combine
  out <- ifelse(x <= x_thresh, poly_part, exp_part)
  
  # Optional: shift everything to have min at 0.5 for higher baseline
  out <- out 
  out <- out +300
  
  
  return(as.vector(scale(out)))
}
xp_prime <- function(x) {
  xprime <- x[-1] - diff(x) / 2
  xpprime <- xprime[-1] - diff(xprime) / 2
  return(xpprime)
}
second_deriv <- function(x, y) {
  xprime <- x[-1] - diff(x) / 2
  dydx <- diff(y) / diff(x)
  d2yd2x <- diff(dydx) / diff(xprime)
  return(d2yd2x)
}
second_deriv_ci <- function(model, data, threshold = "piecewise", score_metric) {
  
  # Predict over a fine sequence of covariate values;
  # avoid the absolute boundaries as often the estimated
  # second derivative may be unstable in those regions
  newcov <- seq(
    quantile(data$covariate, probs = 0.05),
    quantile(data$covariate, probs = 0.95),
    length.out = 500
  )
  newdat <- data.frame(
    covariate = newcov,
    cap = unique(data$cap),
    series = data$series[1]
  )
  
  # We want predictions on the log scale as this is where the
  # function was originally estimated
  preds <- log(
    predict(
      model,
      newdata = newdat,
      type = "link",
      summary = FALSE
    )
  )
  
  # Score the log-function against the true simulated function
  # using an evenly weighted combination of the variogram score and the
  # energy score; this provides a proper multivariate scoring rule that
  # gives better (lower) scores when the predictions capture the true shape
  # of the simulated curve AND the correlations among the true curve's 
  # values along the evaluation points
  preds_zscored <- apply(preds, 1, scale)
  true_curve <- switch(
    threshold,
    piecewise = sim_piecewise(newcov),
    monotonic = sim_monotonic(newcov),
    nonmonotonic = sim_nonmonotonic(newcov),
    cubic = sim_cubic(newcov),
    linear = sim_linear(newcov),
    hyperabundance = sim_poly_then_exp(newcov)
  )
 # score <- (scoringRules::vs_sample(y = true_curve,
 #                                      dat = preds_zscored) +
  #  scoringRules::es_sample(y = true_curve,
  #                              dat = preds_zscored))/2
  if(score_metric == 'vs'){
    score <- scoringRules::vs_sample(y = true_curve, dat = preds_zscored, p=1)
  }
  if(score_metric == 'es'){
    score <- scoringRules::es_sample(y = true_curve, dat = preds_zscored)
  }
  
  # Calculate numerical second derivative
  newcov_prime <- xp_prime(newcov)
  derivs <- matrix(
    NA,
    nrow = NROW(preds),
    ncol = length(newcov_prime)
  )
  for (i in 1:NROW(derivs)) {
    derivs[i, ] <- second_deriv(newcov, preds[i, ])
  }
  
  # Calculate credible intervals of derivatives
  cred <- sapply(
    1:NCOL(derivs),
    function(n) {
      quantile(
        derivs[, n],
        probs = seq(0.1, 0.9, by = 0.1),
        na.rm = TRUE
      )
    }
  )
  
  # Collect data into a data.frame
  deriv_dat <- data.frame(
    newcov_prime,
    lower = cred[1, ],
    med = cred[5, ],
    upper = cred[9, ]
  ) %>%
    # Highlight when credible intervals of 2nd deriv
    # don't overlap with zero
    mutate(
      sig = case_when(
        lower < 0 & upper < 0 ~ 1,
        lower > 0 & upper > 0 ~ 1,
        TRUE ~ 0
      ),
      lower_sig = lower * sig,
      upper_sig = upper * sig
    )
  
  return(list(
    deriv_dat = deriv_dat,
    accuracy = score
  ))
}
threshold_accuracy <- function(model, data, threshold = "piecewise",score_metric) {
  
  # Calculate second derivatives and credible intervals
  deriv_ci <- second_deriv_ci(
    model = model,
    data = data,
    threshold = threshold,
    score_metric = score_metric
  )
  
  return(deriv_ci$accuracy)
}

#### Read Job Array index value into R
args = commandArgs(trailingOnly = TRUE)
slurm = args[1]
slurm = as.numeric(slurm) # imports as character var, not numeric
print(slurm)
# Read in helper parameters 
param_grid <- readRDS("data/param_grid.rds")[slurm,]
design_scenario <- param_grid$design_scenario
simrep <- param_grid$simrep
threshold <- param_grid$threshold

# Read in correct data frame from master list
simdat <- readRDS("data/SimulationInputData.rds")[[design_scenario]][[simrep]]

#### Read external values from SLURM into R
setting = Sys.getenv("SETTING")  # MCMC setting

# MCMC settings, based on assignment above
## Want burn-in to be ~20% of iterations and then thin = (ni - nb) / ideal n.eff (per chain), ideally 30000 in the long one. 
### Assess n.eff via (ni - nb)/nt * nc 
if(setting == "SHORT"){
  ni <- 350;  nt <- 1; nb <- 200; nc <- 4; na = NULL      #quick test to make sure code works
}

form <- switch(param_grid$models,
  'linear' = as.formula(~ covariate),
  'quadratic' = as.formula(~ poly(covariate, 2)),
  'cubic' = as.formula(~ poly(covariate, 3)),
  'GAM' = as.formula(~ s(covariate, k=5)),
  stop("Invalid model type")
)


## specify the number of cores to be uses, should be the same as the model 
options(mc.cores = 4)

##### Run the models!

# Set up the trend_map and fit an N-mixture model
simdat %>%
  # each unique site is a separate process in the SS model
  mutate(trend = as.numeric(site)) %>%
  dplyr::select(trend, series) %>%
  distinct() -> trend_map

# Need some logic to ensure we are running the correct model formulation

# This example uses a thin-plate spline for the
# covariate response
mod<- mvgam(
  formula = 
    obs_count ~ 1,
  trend_formula = form,
  family = nmix(),
  
  # Informative priors for the Intercept (logit(detection probability))
  # and for the Intercept_trend (log(mean latent abundance))
  priors = c(prior(std_normal(),
                   class = Intercept),
             prior(normal(1, 1.5),
                   class = Intercept_trend)),
  data = simdat,
  trend_map = trend_map,
  residuals = FALSE,
  backend = 'rstan',
  
  # parameter testing block
  chains = nc,
  burnin = nb,
  samples = ni,
  thin = nt
)

# Extract key parameters from model output

summary(mod)

es_acc <- threshold_accuracy(model = mod,
                   data = simdat,
                   threshold = threshold,
                   score_metric = 'es')

vs_acc <- threshold_accuracy(model = mod,
                            data = simdat,
                            threshold = threshold,
                            score_metric = 'vs')

# Add this in to our gridline, this is the main output we are looking at
param_grid$vs_accuracy <- vs_acc
param_grid$es_accuracy <- es_acc

saveRDS(param_grid, paste0("results/model_accuracy/model_accuracy", slurm, ".rds"))


