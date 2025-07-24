#!/bin/bash

#SBATCH --job-name=split
#SBATCH --output=split-%j-stdout.log
#SBATCH --error=split-%j-stderr.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --nodes=1

echo "Submitting Slurm job"
singularity exec /idia/software/containers/casa-stable-5.7.0.simg xvfb-run -a casa --nologger --nogui -c split.py
