#### N-Mixture Simulation functions 

# original code by Nick Clark
library(mvgam)
library(gratia)
library(marginaleffects)
library(dplyr)
library(ggplot2); theme_set(theme_classic())
library(patchwork)
library(scoringRules)

# Source custom functions
source('Functions/threshold_functions.R')

# View the types of thresholds available
plot_thresholds()

# Simulate a dataset and plot the true threshold response
threshold <- 'monotonic'
simdat <- sim_nmix(threshold = threshold,
                   n_sites = 25,
                   n_replicates = 4,
                   mean_abundance = 50)

# View the structure of the data
head(simdat, 16)
glimpse(simdat)
simdat <- simdat[[1]][[1]]


ggplot(simdat, aes(x = covariate, y = truth)) +
  geom_point()

 # Set up the trend_map and fit an N-mixture model
simdat %>%
  # each unique site is a separate process in the SS model
  mutate(trend = as.numeric(site)) %>%
  dplyr::select(trend, series) %>%
  distinct() -> trend_map
glimpse(trend_map)

# This example uses a thin-plate spline for the
# covariate response
mod <- mvgam(
  formula = 
    obs_count ~ 1,
  trend_formula = ~
    s(covariate),
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
  backend = 'rstan'
)

# View the model Stan code
stancode(mod)

# Draw the estimated functional response
draw(mod, trend_effects = TRUE)

# Plot estimated latent abundance vs the truth as a function
# of the covariate
plot_predictions(
  mod,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
) +
  geom_point(data = simdat,
             aes(x = covariate,
                 y = truth)) +
  labs(x = 'Covariate',
       y = 'Predicted latent abundance') 

p2<-plot_predictions(
  mod_81,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
) +
  geom_point(data = simdat,
             aes(x = covariate,
                 y = truth)) +
  labs(x = 'Covariate',
       y = 'Predicted latent abundance') +
  geom_line(data = data.frame(
    cov = cov,
    resp = sim_piecewise(x = cov)
  ),
  aes(
    x = cov,
    y = resp
  ))

p1+p2+p3

cov <- seq(-2, 2, length.out = 1000)
p3<-ggplot(
  data.frame(
    cov = cov,
    resp = sim_piecewise(x = cov)
  ),
  aes(
    x = cov,
    y = resp
  )
) +
geom_line(
  linewidth = 1.5,
  col = "darkred")
# Estimated detection probability
avg_predictions(mod, type ='detection')

# Compute "accuracy" of the model wrt whether or not it detected a threshold
# in the correct place. The higher this value (ranges from 0 - 1), the larger
# the overlap between the detected threshold and the true simulated threshold while
# penalising the detection of a threshold outside of the simulated boundary
threshold_accuracy(model = mod,
                   data = simdat,
                   threshold = threshold)

# Plot the estimated 2nd derivative of the covariate response
# function. This plot will highlight any areas where the 2nd derivative
# of the function is strongly supported to be non-zero, indicating possible
# thresholds
plot_function_deriv2(model = mod,
                     data = simdat,
                     threshold = threshold)

# Now try a third order polynomial
mod2 <- mvgam(
  formula = 
    obs_count ~ 1,
  trend_formula = ~
    poly(covariate, 3),
  family = nmix(),
  
  # Informative priors for the Intercept (logit(detection probability)),
  # the Intercept_trend (log(mean latent abundance)) and the 
  # fixed effects
  priors = c(prior(std_normal(),
                   class = Intercept),
             prior(normal(1, 1.5),
                   class = Intercept_trend),
             prior(std_normal(),
                   class = b)),
  data = simdat,
  trend_map = trend_map,
  backend = 'rstan'
)

threshold_accuracy(model = mod2,
                   data = simdat,
                   threshold = threshold)

plot_predictions(
  mod,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
) +
  geom_point(data = simdat,
             aes(x = covariate,
                 y = truth)) +
  labs(x = 'Covariate',
       y = 'Predicted latent abundance')

plot_function_deriv2(model = mod2,
                     data = simdat,
                     threshold = threshold)

# Need to repeat this many times for different types of thresholds and with
# different combinations of n_sites and n_replicates.


second_deriv_ci <- function(model, data, threshold = "piecewise") {
  
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
    cubic = sim_cubic(newcov),
    linear = sim_linear(newcov)
  )
  # small constant to avoid negatives
  score <- (log(scoringRules::vs_sample(y = true_curve,
                                      dat = preds_zscored)+1) +
  log(scoringRules::es_sample(y = true_curve,
                               dat = preds_zscored)+1))/2
  score <- scoringRules::vs_sample(y = true_curve, dat = preds_zscored, p=.5)
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

model = mod_81
data = simdat_81
threshold = 'piecewise'




