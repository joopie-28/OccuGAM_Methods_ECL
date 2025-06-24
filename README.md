# OccuGAMs: Evaluating the performance of GAMs and polynomial approaches to nonlinear occupancy modelling when detection is imperfect

![screenshot](Imagery/HeaderImage.png)

## Last updated: June 24th, 2025

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue?logo=github)](https://github.com/joopie-28/OccuGAM_Methods_ECL)

Methods paper comparing OccuGAMs, AbuGAMs and conventional models for estimation of wildlife occupancy and abundance from cam-trap data.

---

## **Project Background**

This GitHub repository includes all code to implement the analyses featured in Sassen et al.'s OccuGAMs analysis, which comprises a comparison of polynomial and penalised spline (i.e. GAM) approaches to modelling non-linear responses to disturbance covariates. The repository also includes additional annotated code templates to implement 'OccuGAMs' in both STAN and JAGS.

The repository contains 2 types of scripts:
-**R scripts** consitute the majority of scripts. They are used to i) process data into formats suitable for modelling, ii) create simulated datasets, iii) Fit models to data - we fit 1000s of complex models utilising a High Performance Computing cluster [(HPC)](https://rcc.uq.edu.au/systems/high-performance-computing/bunya) and iv) Gather processed models from HPC cluster and visualise the results.
-**SLURM scripts** are used for batch processing of analyses on the HPC environment.




