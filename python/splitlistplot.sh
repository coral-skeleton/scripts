#!/bin/bash

#SBATCH --job-name='splitlistplot'
#SBATCH --output=%x-%j-stdout.log
#SBATCH --error=%x-%j-stderr.log
#SBATCH --mail-user=jc@saao.ac.za
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=03:00:00
#SBATCH --nodes=1

echo "Submitting Slurm job"
singularity exec /idia/software/containers/casa-stable-5.7.0.simg xvfb-run -a casa --nologger --nogui -c splitlistplot.py