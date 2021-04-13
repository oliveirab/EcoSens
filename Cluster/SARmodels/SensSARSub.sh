#!/bin/bash
#SBATCH --time=24:00:00


for i in {1..12} 
do

	JOB_NAME="M${i}_SensSAR"
	ER_NAME="/home/brolivei/SensProd/SARmodels/slurm-log/M${i}-%j-err.txt"
	OUT_NAME="/home/brolivei/SensProd/SARmodels/slurm-log/M${i}-%j-out.txt"

	sbatch --job-name=$JOB_NAME -e $ER_NAME -o $OUT_NAME SensSARJob.sh $i


done