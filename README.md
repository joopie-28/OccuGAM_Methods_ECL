# OccuGAMs: non-linear occupancy and abundance modelling with imperfect detection

![screenshot](Imagery/HeaderImage.png)

## Last updated: August 22nd, 2025

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue?logo=github)](https://github.com/joopie-28/OccuGAM_Methods_ECL)

Methods paper comparing OccuGAMs, AbuGAMs and conventional models for estimation of wildlife occupancy and abundance from cam-trap data.

---

## **Project Background**

This GitHub repository includes all code to implement the analyses featured in Sassen et al.'s OccuGAMs analysis, which comprises a comparison of polynomial and penalised spline (i.e. GAM) approaches to modelling non-linear responses to disturbance covariates. The repository also includes additional annotated code templates to implement 'OccuGAMs' in both STAN and JAGS.

The repository contains 2 types of scripts:
- **R scripts** consitute the majority of scripts. They are used to i) process data into formats suitable for modelling, ii) create simulated datasets, iii) Fit models to data - we fit 1000s of complex models utilising a High Performance Computing cluster [(HPC)](https://rcc.uq.edu.au/systems/high-performance-computing/bunya) and iv) Gather processed models from HPC cluster and visualise the results.
- **SLURM scripts** are used for batch processing of analyses on the HPC environment.

The camera trap data used in this analysis was prepared using a multi-step data cleaning pipeline that is backed up on GitHub but is currently private due to data-sharing agreements with our collaborators. The dataset provided in this repository has been de-identified, with all latitude and longitude coordinates removed to protect sampling locations.


## **Repository Structure**

The repository is organised into 5 main folders, of which the folders **Analyses**, **Functions**, **Inputs** are the 3 folders strictly required to reproduce all analyses. The folder **Outputs** contains the majority of results including plots and model summaries. **ModelTemplates** contains example code in JAGS and STAN for integrating penalised splines in occupancy and n-mixture models.

### 0. Functions and Pre-Requisites
- Functions are stored separately in this repository, and are called at the start of each separate script to promote readability and modularity of code. The only exception is in the HPC Packages, as these require all code to be bundled.
- All functions are stored in the `Functions/` subfolder.

### 1. Data preparation
- `Analyses/R/01_Analysis_CreateDataBundles.R` contains all code to generate count and detection history matrices from ECL camera trap data. The outputs are in the 'unmarkedFrame' format which conveniently stores all relevant detections and covariates. These processed files are stored in the `Outputs/UMF.List` subfolder. A copy of these umf.list is also stored in the HPC packages (see Model Fitting on the HPC environment).
- `Analyses/Simulations/02_Simulations_BuildSimFramework.R` contains all code to simulate 2,400 datasets for the simulation portion of the study. We note that there is a stochastic component to this and thus newly simulated datasets will deviate from the ones produced in this work. The simulated dataset used to fit models for this paper is included in the HPC packages (see Model Fitting on the HPC environment).

### 2. Model Fitting on the HPC environment
- `Analyses/HPC_Packages` contains completed code and data packages that were sent/can be sent to a HPC environment. 
- `Analyses/HPC_Packages/05_HPC_Comprehensive` contains all code, data for the case study portion of the paper. This comprises fitting linear, quadratic, cubic and GAM occupancy and n-mixture models for 4 species of tropical mammals (*Sus scrofa*, *Rusa unicolor*, *Macaca nemestrina*, genus *Muntiacus*).
- `Analyses/HPC_Packages/06_HPC_Simulations` contains all code, data for the simulation portion of the paper. This comprises fitting linear, quadratic, cubic and GAM occupancy and n-mixture models to 2,400 different scenario datasets.

### 3. Analysing and visualising results
Visualising of results was performed on a local device and required extracting the key parameters and information from the HPC_packages. 
- `Analyses/R/03_Process_Results.R` contains all code to produce the tables and graphs from the case study portion of the manuscript. The resulting artefacts are stored in `Outputs/`. 
- `Analyses/Simulations/03_Process_Simulation_Results.R` contains all code to produce the tables and graphs from the simulation portion of the manuscript. The resulting artefacts are stored in `Outputs/`. 

### 4. Final results and figures
The `Outputs/` subfolder contains all study results. We note that the raw, unprocessed results (i.e. raw models) were too large to host on GitHub and hence are provided in the accompanying FigShare repository: [Figshare Results Repository](https://figshare.com/projects/OccuGAMs_non-linear_occupancy_and_abundance_modelling_with_imperfect_detection/254597).
- `Outputs/Abundance_Plots` contains the habitat association graphs (Fig. 3, Fig S1-S3 in main report). These are modular and were stitched together as a single graph in powerpoint. It includes both occupancy and abundance plots.
- `Outputs/LOO Performance Plots` contains the predictive performance graphs as measured by the ELPD (Fig. 4, Fig S5).
- `Outputs/Main Figures` contains the NRMSE graphs (Fig. 5, Fig S4).
- `Outputs/Main Results` contains the model PPE tables (Table. S1-S2).
- `Outputs/Simulations` contains the simulation study graphs (Fig. 1-2). These are modular and were stitched together as single graphs in powerpoint.
- `Outputs/UMF.List` contains the unmarkedFrame data files (see Data Preparation)

## **OccuGAM model code in JAGS and STAN**
- In order to help academic and applied ecologists alike, we provide clear model code extracts to integrate penalised splines into both occupancy and n-mixture models in both JAGS and STAN. 
- We are only aware of prior use of this approach in occupancy models (i.e. not abundance/n-mixture) and JAGS (i.e. not STAN)
- Our hope is that this will facilitate greater model uptake, particularly in (relative) abundance estimation, given the benefits compared to polynomials.
- `ModelCodeTemplates/JAGS` contains the code to run occupancy and n-mixture models in JAGS. We note that this code is dependent on packages such as `mgcv` and for initial set-up. We refer the user to the `Analyses/HPC_Packages/05_HPC_Comprehensive/Abundance/code/01_HPC_FitPolynomialAbuJags.R` file if they wish to understand this further.
- `ModelCodeTemplates/STAN` contains the code to run n-mixture models in STAN. We note that this code is dependent on packages such as `mgcv` and for initial set-up. We also point the user to the `mvgam` r package, which we used to construct the STAN n-mixture models, and the `Flocker` r package for occupancy modelling. They offer an intuitive formula-based interface to STAN and provide an easy entry point.

