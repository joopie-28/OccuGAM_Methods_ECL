# Simulate a piecewise threshold response to a standardized covariate
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

# Simulate a monotonic response to a standardized covariate
sim_monotonic <- function(x) {
  out <- exp(2.2 * x) / (6 + exp(2 * x)) * -1
  return(as.vector(scale(out)))
}

# Simulate a mnon-onotonic response to a standardized covariate
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

# Simulate a cubic response to a standardized covariate
sim_cubic <- function(x) {
  out <- -0.5 * x^2 + 1.5 * x^3 * -1
  return(as.vector(scale(out)))
}

# Simulate a linear response to a standardized covariate
sim_linear <- function(x) {
  out <- 0.4 * x
  return(as.vector(scale(out)))
}

# hyperabundance
sim_poly_then_exp <- function(x) {
  x_thresh <- quantile(x, probs = 0.5)
  
  # Centered covariate for smooth polynomial dip
  poly_part <- 3.2 * (x - x_thresh)^4  # Quartic, starts at ~0.5, dips
  
  # Exponential increase after threshold
  exp_part <- exp(3 * (x - x_thresh)) - 1
  exp_vals <- exp_part[x > x_thresh]
  
  # Combine
  out <- ifelse(x <= x_thresh, poly_part, exp_part)
  
  # Optional: shift everything to have min at 0.5 for higher baseline
  out <- out - min(out)
  out <- out + 50
  
  return(as.vector(scale(out)))
}

# Some boundaries of where thresholds are, just for plotting
# purposes
threshold_boundaries <- function(threshold = "piecewise") {
  min_bound <- switch(threshold,
                      piecewise = -0.75,
                      monotonic = 0,
                      cubic = -1.5,
                      linear = Inf
  )
  
  max_bound <- switch(threshold,
                      piecewise = -0.35,
                      monotonic = 1,
                      cubic = -0.5,
                      linear = Inf
  )
  
  min_bound2 <- switch(threshold,
                       piecewise = Inf,
                       monotonic = Inf,
                       cubic = 0.8,
                       linear = Inf
  )
  
  max_bound2 <- switch(threshold,
                       piecewise = Inf,
                       monotonic = Inf,
                       cubic = 1.3,
                       linear = Inf
  )
  
  return(list(
    min_bound = min_bound,
    max_bound = max_bound,
    min_bound2 = min_bound2,
    max_bound2 = max_bound2
  ))
}

# Function to plot the different types of threshold responses available
plot_thresholds <- function() {
  cov <- seq(-2, 2, length.out = 1000)
  wrap_plots(
    ggplot(
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
        col = "darkred"
      ) +
      labs(
        x = "Response",
        y = "Covariate",
        title = "Piecewise"
      ),
    ggplot(
      data.frame(
        cov = cov,
        resp = sim_nonmonotonic(x = cov)
      ),
      aes(
        x = cov,
        y = resp
      )
    ) +
      geom_line(
        linewidth = 1.5,
        col = "darkred"
      ) +
      labs(
        x = "Response",
        y = "Covariate",
        title = "Nonmonotonic"
      ),
    ggplot(
      data.frame(
        cov = cov,
        resp = sim_monotonic(x = cov)
      ),
      aes(
        x = cov,
        y = resp
      )
    ) +
      geom_line(
        linewidth = 1.5,
        col = "darkred"
      ) +
      labs(
        x = "Response",
        y = "Covariate",
        title = "Monotonic"
      ),
    ggplot(
      data.frame(
        cov = cov,
        resp = sim_cubic(x = cov)
      ),
      aes(
        x = cov,
        y = resp
      )
    ) +
      geom_line(
        linewidth = 1.5,
        col = "darkred"
      ) +
      labs(
        x = "Response",
        y = "Covariate",
        title = "Cubic"
      ),
    ggplot(
      data.frame(
        cov = cov,
        resp = sim_linear(x = cov)
      ),
      aes(
        x = cov,
        y = resp
      )
    ) +
      geom_line(
        linewidth = 1.5,
        col = "darkred"
      ) +
      labs(
        x = "Response",
        y = "Covariate",
        title = "Linear"
      ),
    ggplot(
      data.frame(
        cov = cov,
        resp = sim_poly_then_exp(x = cov)
      ),
      aes(
        x = cov,
        y = resp
      )
    ) +
      geom_line(
        linewidth = 1.5,
        col = "darkred"
      ) +
      labs(
        x = "Response",
        y = "Covariate",
        title = "Linear"
      ),
    ncol = 2
  )
}

# Function to validate that a numeric argument is a proportion
validate_proportional <- function(x) {
  s <- substitute(x)
  x <- base::suppressWarnings(as.numeric(x))
  if (length(x) != 1L || anyNA(x)) {
    stop("Argument '", s, "' must be a single numeric value", call. = FALSE)
  }

  if (x < 0 || x > 1) {
    stop(
      "Argument '",
      s,
      "' must be a proportion ranging from 0 to 1, inclusive",
      call. = FALSE
    )
  }
}

# Function to validate that an argument is a positive integer
validate_pos_integer <- function(x) {
  s <- substitute(x)
  x <- base::suppressWarnings(as.numeric(x))
  if (length(x) != 1L || anyNA(x)) {
    stop("Argument '", s, "' must be a single numeric value", call. = FALSE)
  }

  if (sign(x) != 1) {
    stop("Argument '", s, "' must be a positive integer", call. = FALSE)
  } else {
    if (x %% 1 != 0) {
      stop("Argument '", s, "' must be a positive integer", call. = FALSE)
    }
  }
}

# Function to simulate imperfect observations of a true abundance
# that varies in response to a single covariate
sim_nmix <- function(
    threshold = c(
      "piecewise",
      "monotonic",
      'nonmonotonic',
      "cubic",
      "linear",
      "hyperabundance"
    ),
    p_detect = 0.5,
    mean_abundance = 20,
    n_sites = 25,
    n_replicates = 4) {
  
  # Validate arguments
  threshold <- match.arg(threshold)
  validate_proportional(p_detect)
  validate_pos_integer(n_sites)
  validate_pos_integer(n_replicates)
  validate_pos_integer(mean_abundance)
  if (mean_abundance < 8) {
    stop(
      "Not recommended to fit N-mixtures when mean abundance < 8",
      call. = FALSE
    )
  }

  # Simulate the environmental covariate
  sim_cov <- runif(n = n_sites, min = -2, max = 2)

  # Simulate the abundance response over the n_sites
  cov_response <- switch(
    threshold,
    piecewise = sim_piecewise(sim_cov),
    monotonic = sim_monotonic(sim_cov),
    nonmonotonic = sim_nonmonotonic(sim_cov) ,
    cubic = sim_cubic(sim_cov),
    linear = sim_linear(sim_cov),
    hyperabundance = sim_poly_then_exp(sim_cov)
  )

  # Simulate the latent abundances
  latent_n <- rpois(
    n = n_sites,
    lambda = pmax(
      1,
      mean_abundance +
        cov_response * (mean_abundance / 2.5) + rnorm(n=n_sites, mean = 0, sd = (mean_abundance*0.2))
    )
  )

  # Take n_replicates imperfect observations at each site
  obs_dat <- do.call(rbind, lapply(1:n_sites, function(x) {
    data.frame(
      site = x,
      replicate = 1:n_replicates,
      obs_count = rbinom(
        n = n_replicates,
        size = latent_n[x],
        prob = p_detect
      ),
      covariate = sim_cov[x],
      truth = latent_n[x],
      threshold = threshold,
      p_detect = p_detect
    )
  })) %>%
    
    # Finish wrangling for mvgam
    mutate(
      site = factor(site),
      replicate = factor(replicate),
      cap = max(latent_n) * 2
    ) %>%
    mutate(series = factor(paste0(
      "site_", site,
      "_rep_", replicate
    ))) %>%
    dplyr::select(
      series, site:covariate, cap,
      truth, threshold, p_detect
    )

  # Return the simulated data
  return(obs_dat)
}

# Functions to calculate a second derivative numerically
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

# Function to compute second derivatives and score whether or not
# a threshold was detected appropriately
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
  score <- log(scoringRules::vs_sample(y = true_curve,
                                       dat = preds_zscored)) *
    log(scoringRules::es_sample(y = true_curve,
                                dat = preds_zscored))

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

threshold_accuracy <- function(model, data, threshold = "piecewise") {
  
  # Calculate second derivatives and credible intervals
  deriv_ci <- second_deriv_ci(
    model = model,
    data = data,
    threshold = threshold
  )
  
  return(deriv_ci$accuracy)
}

# Function to numerically compute the 2nd derivative of a
# covariate response function, and return a ggplot of this
plot_function_deriv2 <- function(model, data, threshold = "piecewise") {
  
  # Calculate second derivatives and credible intervals
  deriv_ci <- second_deriv_ci(
    model = model,
    data = data,
    threshold = threshold
  )
  plot_dat <- deriv_ci$deriv_dat

  # Plot
  p <- ggplot(
    plot_dat,
    aes(
      x = newcov_prime,
      y = med
    )
  ) +
    geom_ribbon(
      aes(
        ymax = upper,
        ymin = lower
      ),
      alpha = 0.2
    ) +
    geom_ribbon(
      aes(
        ymax = upper_sig,
        ymin = lower_sig
      ),
      fill = "#B97C7C"
    ) +
    geom_line(linewidth = 0.75) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      y = "f''(covariate)",
      x = "covariate",
      title = paste0(threshold, " threshold score = ",
                     round(deriv_ci$accuracy, 4))
    )

  # Add boundaries for the true threshold
  boundaries <- threshold_boundaries(threshold = threshold)
    p <- p +
      geom_vline(
        xintercept = boundaries$min_bound,
        linetype = "dashed"
      ) +
      geom_vline(
        xintercept = boundaries$min_bound2,
        linetype = "dashed"
      ) +
      geom_vline(
        xintercept = boundaries$max_bound,
        linetype = "dashed"
      ) +
      geom_vline(
        xintercept = boundaries$max_bound2,
        linetype = "dashed"
      )

  return(p)
}
