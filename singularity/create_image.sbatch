#!/bin/bash
#SBATCH --job-name=build_icefinder
#SBATCH --output=build_icefinder.log
#SBATCH --error=build_icefinder.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH -p interactive

ml apptainer
apptainer build icefinder.sif singularity.def
