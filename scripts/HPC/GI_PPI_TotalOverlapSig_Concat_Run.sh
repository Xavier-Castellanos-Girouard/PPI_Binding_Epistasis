#!/bin/bash

#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/xaviercg/scratch/log/%j.out


source ~/scratch/Yeast_Epistasis_Energetics_venv/bin/activate

python ~/scratch/scripts/HPC_Module_Overlap_Analysis_Concat.py

deactivate
