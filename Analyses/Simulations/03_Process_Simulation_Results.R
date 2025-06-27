#### Script for processing the 9,600 mvgam simulations
 
library(purrr)
library(ggplot2)
library(dplyr)
library(mvgam)
library(tidybayes)
library(patchwork)

# import the data 

rds.files <- list.files("/Users/sassen/Desktop/06_HPC_Simulations/results/model_accuracy")

# combine into single frame
combined_df <- map_dfr(paste0("/Users/sassen/Desktop/06_HPC_Simulations/results/model_accuracy/",
                              rds.files), 
                       readRDS)

# build the visuals

# Define the desired order for 'threshold'
threshold_order <- c("linear", "cubic", "monotonic", "piecewise", 'hyperabundance')

# Filter and reorder data
filtered_df <- combined_df %>%
  filter(scenario == 'large') |>
  mutate(threshold = factor(threshold, levels = threshold_order))

# Create the 16-panel figure
p<-ggplot(filtered_df, aes(x = models, y = vs_accuracy)) +
  geom_boxplot() +
  facet_grid(rows = vars(threshold), cols = vars(species)) +
  theme_minimal() +
  labs(x = "Models", y = "Accuracy", title = "Model Accuracy by Species and Threshold") +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))  # Rotate x-axis labels

# Print the plot
print(p)

pdf("/Users/sassen/Desktop/SimResults.pdf", width=10, height = 10)
print(p)
dev.off()

# Rank plots


# Assuming your data is in a data frame called 'df'
vs_result <- combined_df |> 
  filter(scenario == 'large') |>
#  filter(threshold != 'linear') |>
  group_by(design_scenario, simrep) |>
  arrange(vs_accuracy, .by_group = TRUE)|>
  mutate(rank = row_number()) |>
  ungroup() |>
  mutate(models = factor(models, levels = c('linear', 'quadratic', 'cubic', 'GAM')))
  

es_result <- combined_df |> 
  filter(scenario == 'large') |>
 # filter(threshold != 'linear') |>
  group_by(design_scenario, simrep) |>
  arrange(es_accuracy, .by_group = TRUE)|>
  mutate(rank = row_number()) |>
  ungroup() |>
  mutate(models = factor(models, levels = c('linear', 'quadratic', 'cubic', 'GAM')))


pdf("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Simulations/SimResults_vs.pdf", width=10, height = 5)

ggplot(data = vs_result, aes(x = rank, y = models, fill = models)) +
  geom_violin(bw = 0.25) +
  facet_wrap(~threshold) +
  #geom_jitter(width = 0.1, height = 0.1, alpha = 0.1, size = 0.1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    width = 0.2,
    color = "black",
    fatten = 2
  ) + ylab('Model Formulation') + xlab('Variogram Rank')+
  scale_fill_manual(values = c(
    "linear" = "lightgreen",
    "quadratic" = 'steelblue',
    "cubic" = "purple",
    "GAM" = "darkgrey"
  ),
  drop = FALSE
  ) + labs(fill = "Model")+
  
  theme_classic() +theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
                         strip.background = element_rect(color = "black", size = 1), 
                         strip.text = element_text(size = 16),axis.title = element_text(size=14), 
                         axis.text = element_text(size=12))# Rotate x-axis labels
dev.off()


pdf("/Users/sassen/OccuGAM_Methods_ECL/Outputs/Simulations/SimResults_es.pdf", width=10, height = 5)

ggplot(data = es_result, aes(x = rank, y = models, fill = models)) +
  geom_violin(bw = 0.25) +
  facet_wrap(~threshold) +
  #geom_jitter(width = 0.1, height = 0.1, alpha = 0.1, size = 0.1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    width = 0.2,
    color = "black",
    fatten = 2
  ) + ylab('Model Formulation') + xlab('Energy Rank')+
  scale_fill_manual(values = c(
    "linear" = "lightgreen",
    "quadratic" = 'steelblue',
    "cubic" = "purple",
    "GAM" = "darkgrey"
  ),
  drop = FALSE
  ) + labs(fill = "Model")+
  theme_classic() +theme(panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
                         strip.background = element_rect(color = "black", size = 1), 
                         strip.text = element_text(size = 16),axis.title = element_text(size=14), 
                         axis.text = element_text(size=12))# Rotate x-axis labels

dev.off()

#### Create Threshold plot to show how we did the simulations

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
        x = "",
        y = "",
        title = ""
      ) + 
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            strip.background = element_rect(color = "black", size = 1), 
            strip.text = element_text(size = 16),axis.title = element_text(size=14), 
            axis.text = element_text(size=12)),# Rotate x-axis labels,
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
        x = "",
        y = "",
        title = ""
      ) + theme_void() ,
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
plot_thresholds()


pdf('Outputs/Simulations/SimulationCurves.pdf', height= 10, width = 15)
cov <- seq(-2, 2, length.out = 1000)
p <- ggplot(
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
    x = "",
    y = "",
    title = ""
  ) + theme_void()

print(p)

p <- ggplot(
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
    x = "",
    y = "",
    title = ""
  ) + theme_void()

print(p)


p <- ggplot(
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
    x = "",
    y = "",
    title = ""
  ) + theme_void()

print(p)

p <- ggplot(
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
    x = "",
    y = "",
    title = ""
  ) + theme_void()

print(p)

p <- ggplot(
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
    x = "",
    y = "",
    title = ""
  ) + theme_void()

print(p)

p <- ggplot(
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
    x = "",
    y = "",
    title = ""
  ) + theme_void()

print(p)

dev.off()

#### Create an example plot from our real simulated datasets

simdat <- readRDS("/Users/sassen/Dropbox/Joop MSc QBIO project/Joop OccuGAMs methods paper/OccuGAMs Methods - Results/Simulation Study/06_HPC_Simulations/data/SimulationInputData.rds") 

# Select a random dataset from the monotonic scenarios
simdat.ex <- simdat[[9]][[5]]

pdf('Outputs/Simulations/ExampleSimulationSet.pdf', height= 3, width = 4.7)
# Plot the simulated aata
ggplot(simdat.ex, aes(x = covariate, y = truth)) +
  geom_point(shape = 21, fill = "darkred", color = "black", stroke = 0.3) + theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),       # plot area stays white
        plot.background = element_rect(fill = NA, color = NA), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        axis.line = element_line(linewidth = .0),
                                     strip.background = element_rect(color = "black", size = .5), 
                                     strip.text = element_text(size = 16, color = 'black'),axis.title = element_text(size=14), 
                                     axis.text = element_text(size=12, color = 'black')) +
  labs(
    x = 'Covariate',
    y = 'Latent N'
  )

dev.off()

## some key stats

vs_result |>
  group_by(threshold, models) |>
  summarise(Mean_ES_RANK = mean(rank),
            Mean_ES = mean(vs_accuracy)) |> View()


#### Create additional example with fitted curves for a hyperabundance scenario

simdat.hyper <- simdat[[10]][[16]]

# Plot the simulated aata
ggplot(simdat.hyper, aes(x = covariate, y = truth)) +
  geom_point(shape = 21, fill = "darkred", color = "black", stroke = 0.3) + theme_classic()+
  theme(panel.background = element_rect(fill = "white", color = NA),       # plot area stays white
        plot.background = element_rect(fill = NA, color = NA), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        axis.line = element_line(linewidth = .0),
        strip.background = element_rect(color = "black", size = .5), 
        strip.text = element_text(size = 16, color = 'black'),axis.title = element_text(size=14), 
        axis.text = element_text(size=12, color = 'black')) +
  labs(
    x = 'Covariate',
    y = 'Latent N'
  )

# Set up the trend_map and fit an N-mixture model
simdat.hyper %>%
  # each unique site is a separate process in the SS model
  mutate(trend = as.numeric(site)) %>%
  dplyr::select(trend, series) %>%
  distinct() -> trend_map

# Set parameters
ni <- 350;  nt <- 1; nb <- 200; nc <- 4; na = NULL  

options(mc.cores = 4)

full.mod.list<- lapply(list('linear', 'quadratic', 'cubic', 'GAM'), FUN = function(mod){
  # Pick the correct model
  form <- switch(mod,
                 'linear' = as.formula(~ covariate),
                 'quadratic' = as.formula(~ poly(covariate, 2)),
                 'cubic' = as.formula(~ poly(covariate, 3)),
                 'GAM' = as.formula(~ s(covariate, k=5)),
                 stop("Invalid model type")
  )
  
  mod <- mvgam(
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
    data = simdat.hyper,
    trend_map = trend_map,
    residuals = FALSE,
    backend = 'rstan',
    
    # parameter testing block
    chains = nc,
    burnin = nb,
    samples = ni,
    thin = nt
  )
  return(mod)
})

names(full.mod.list) <- c('linear', 'quadratic', 'cubic', 'GAM')

## Plot the examples
curves.linear <- mvgam::predictions(
  full.mod.list$linear,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
)

curves.quad <- mvgam::predictions(
  full.mod.list$quadratic,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
)

curves.cubic <- mvgam::predictions(
  full.mod.list$cubic,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
)

curves.gam <- mvgam::predictions(
  full.mod.list$GAM,
  newdata = datagrid(covariate = seq(-2, 2, length.out = 100)),
  by = "covariate",
  type = "link"
)

pdf('Outputs/Simulations/ExampleSimulationNonMono_G.pdf', height= 5, width = 7)
p <-  ggplot() + 
  geom_line(data=dftest,
            aes(x = cov, y=resp),
            colour="darkred") +
    geom_line(data = curves.linear, 
              aes(x = covariate, y=estimate), 
              colour='green') + 
    geom_ribbon(data= curves.linear, aes(x = covariate, ymin = conf.low, ymax = conf.high),
               alpha =0.2, 
                fill = 'green')+
  geom_line(data = curves.quad, 
            aes(x = covariate, y=estimate), 
            colour='steelblue') + 
  geom_ribbon(data= curves.quad, aes(x = covariate, ymin = conf.low, ymax = conf.high),
             alpha =0.2, 
              fill = 'steelblue')+
  geom_line(data = curves.cubic, 
            aes(x = covariate, y=estimate), 
            colour='purple') + 
  geom_ribbon(data= curves.cubic, aes(x = covariate, ymin = conf.low, ymax = conf.high),
            alpha =0.2, 
             fill = 'purple')+
  geom_line(data = curves.gam, 
            aes(x = covariate, y=estimate), 
            colour='black') + 
  geom_ribbon(data= curves.gam, aes(x = covariate, ymin = conf.low, ymax = conf.high),
              alpha =0.2, 
              fill = 'black') +
  geom_point(data = simdat.hyper,
               aes(x = covariate,
                   y = truth),shape = 21, fill = "black",stroke = 0.5) + 
    geom_point(data = simdat.hyper,
              aes(x = covariate,
                  y = obs_count),
              shape = 21, fill = "darkred",stroke = 0.5, size=0.2) +
    labs(x = 'Covariate',
         y = 'Abundance') +
    theme(panel.background = element_rect(fill = "white", color = NA),       # plot area stays white
          plot.background = element_rect(fill = NA, color = NA), 
          panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
          axis.line = element_line(linewidth = .0),
          strip.background = element_rect(color = "black", size = .5), 
          strip.text = element_text(size = 16, color = 'black'),axis.title = element_text(size=14), 
          axis.text = element_text(size=12, color = 'black'))
print(p)
dev.off()

pdf('Outputs/Simulations/NonMonoCurve.pdf', height = 5, width = 7)

cov <- seq(-2, 2, length.out = 1000)
p <- ggplot(
  data.frame(
    cov = cov,
    resp = sim_nonmonotonic(x = cov)*(100/2.5)+100
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
    x = "",
    y = "",
    title = ""
  ) 

print(p)

dev.off()

