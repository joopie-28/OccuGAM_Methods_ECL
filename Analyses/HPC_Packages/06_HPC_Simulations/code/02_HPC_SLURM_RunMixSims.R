#!/bin/bash
# SPECIFY YOUR GROUP NAME
#SBATCH --account=a_luskin_ecl

## SPECIFY YOUR COMPUTATIONAL REQUIREMENTS 

# Select 1 node per job
#SBATCH --nodes=1

# Select 1 task per node (b/c R is not MPI)
#SBATCH --ntasks=1

# Select 4 CPUS per node (one per MCMC chain)
#SBATCH --cpus-per-task=4

# Select 50000 MB (50 GB) of memory per node 
#SBATCH --mem=50000

# Ensure we are in the general queue, not AI, debug, or GPU
#SBATCH --partition=general

# Select 4 hours (h:m:s format) of walltime for small mods; we should not need that long
#SBATCH --time=4:00:00

# SPECIFY THE JOB ARRAY-
#SBATCH --array=1-9600

# SPECIFY THE JOB NAME
#SBATCH --job-name=NmixSims

# SPECIFY .err AND .out FILE LOCATIONS
#SBATCH --output=/scratch/user/s4633921/06_HPC_Simulations/results/OUT/slurm-%A_%a.out
#SBATCH --error=/scratch/user/s4633921/06_HPC_Simulations/results/ERRORS/slurm-%A_%a.err

# SET THE 'setting' VARIABLE THAT WILL BE LOADED IN R
export SETTING="SHORT" 

# SPECIFY THE PBS WORKING DIRECTORY AND PRINT TO VERIFY
cd $SLURM_SUBMIT_DIR
pwd

# LOAD THE R SCRIPT, SUBMIT JOB AND SPECIFY THE ARRAY INDEX 
apptainer exec /scratch/user/s4633921/stan_simulation_container.sif Rscript /scratch/user/s4633921/06_HPC_Simulations/code/01_RunSimulationModels.R $SLURM_ARRAY_TASK_ID





