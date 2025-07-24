#!/bin/bash

#SBATCH --job-name='caracal'
#SBATCH --output=%x-%j-stdout.log
#SBATCH --error=%x-%j-stderr.log
#SBATCH --mail-user=jc@saao.ac.za
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=232GB
#SBATCH --time=150:00:00
#SBATCH --nodes=1

echo "Submitting Slurm job"
source /scratch3/projects/meerchoirs/jcviljoen/caracal-new/caracalenv/bin/activate
caracal -c full.yml -ct singularity -sid /software/astro/caracal/STIMELA_IMAGES_1.7.7
deactivate

# https://github.com/caracal-pipeline/caracal/issues/625n