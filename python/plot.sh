#!/bin/bash

#SBATCH --job-name=plot
#SBATCH --output=plot-%j-stdout.log
#SBATCH --error=plot-%j-stderr.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=10:00:00
#SBATCH --nodes=1

echo "Submitting Slurm job"
singularity exec /idia/software/containers/casa-stable-5.7.0.simg xvfb-run -a casa --nologger --nogui -c plot.py
