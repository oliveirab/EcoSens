#!/bin/bash
#SBATCH --time=3-24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --partition=bmm
#SBATCH --job-name=Imputation
#SBATCH -D /home/brolivei/SensProd/Imputation/
#SBATCH -e /home/brolivei/SensProd/Imputation/stderr-%j.txt
#SBATCH -o /home/brolivei/SensProd/Imputation/stdout-%j.txt
#SBATCH --mail-user=brolivei@ucdavis.edu
#SBATCH --mail-type=ALL

module load R
Rscript imputationCode.R Rcodeoutput.txt