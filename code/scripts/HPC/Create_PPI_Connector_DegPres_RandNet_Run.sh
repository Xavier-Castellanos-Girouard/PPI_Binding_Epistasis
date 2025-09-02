#!/bin/bash

#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=/home/xaviercg/scratch/log/%j.out



source ~/scratch/Yeast_Epistasis_Energetics_venv/bin/activate

python ~/scratch/scripts/HPC_Create_PPI_Connector_DegPres_RandNet.py


deactivate
