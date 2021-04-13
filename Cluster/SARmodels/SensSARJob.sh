#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1 # jobs running in a single node
#SBATCH --mem=0
#SBATCH -D /home/brolivei/SensProd/SARmodels/
#SBATCH --mail-user=brolivei@ucdavis.edu
#SBATCH --mail-type=ALL

i=$1

module load R
Rscript SensSARCode.R $i